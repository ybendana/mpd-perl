package MPD::Bed;

# ABSTRACT: This package manipulates Bed Files

use 5.10.0;

use Moose 2;
use MooseX::Types::Path::Tiny qw/ AbsPath AbsFile File /;
use namespace::autoclean;

use Type::Params qw/ compile /;
use Types::Standard qw/ :types /;
use Scalar::Util qw/ reftype /;
use Path::Tiny;
use Try::Tiny;

use Data::Dump qw/ dump /; # for debugging

use MPD::Bed::Raw;
use DDP;

use Scalar::Util qw/looks_like_number/;

with 'MPD::Role::Message';

our $VERSION = '0.001';

has BedFile => ( is => 'ro', isa => AbsPath, coerce => 1 );
has Entries => (
  traits  => ['Array'],
  is      => 'ro',
  isa     => 'ArrayRef[MPD::Bed::Raw]',
  handles => {
    all_entries   => 'elements',
    no_entires    => 'is_empty',
    count_entries => 'count',
  },
  default => sub { [] },
);

has CoveredChr => (
  traits  => ['Hash'],
  is      => 'ro',
  isa     => 'HashRef',
  handles => {
    exists_chr => 'exists',
    all_chrs   => 'keys',
    get_chr    => 'get',
    no_chrs    => 'is_empty',
    count_chrs => 'count',
  },
  default => sub { {} },
);

has CoveredSite => (
  traits  => ['Hash'],
  is      => 'ro',
  isa     => 'HashRef',
  handles => {
    exists_site => 'exists',
    all_sites   => 'keys',
    get_site    => 'get',
    no_sites    => 'is_empty',
    count_sites => 'count',
  },
  default => sub { {} },
);

sub Entries_as_BedFileLetter {
  my $self = shift;

  my @strs;

  my $entriesAref = $self->Entries_as_aref();

  for my $e (@$entriesAref) {
    my $line;
    my $chr = shift @$e;
    if ( $chr == 23 ) {
      $line = "chr" . join( "\t", ( 'M', @$e ) );
    }
    elsif ( $chr == 24 ) {
      $line = "chr" . join( "\t", ( 'X', @$e ) );
    }
    elsif ( $chr == 25 ) {
      $line = "chr" . join( "\t", ( 'Y', @$e ) );
    }
    else {
      $line = "chr" . join( "\t", ( $chr, @$e ) );
    }
    push @strs, $line;
  }
  return join( "\n", @strs );
}

sub Entries_as_BedFile {
  my $self = shift;

  my @strs;

  my $entriesAref = $self->Entries_as_aref();
  for my $e (@$entriesAref) {
    my $line = "chr" . join( "\t", @$e );
    push @strs, $line;
  }
  return join( "\n", @strs );
}

sub Entries_as_aref {
  my $self = shift;

  my @array;

  for my $e ( $self->all_entries ) {
    push @array, $e->as_aref;
  }
  return \@array;
}

# SiteNames gives a hash ref of all sites in the bed object and their
# corresponding names
sub SiteNames {
  my $self = shift;

  my %hash;

  for my $e ( $self->all_entries ) {
    my $chr = $e->Chr;
    for ( my $i = $e->Start; $i <= $e->End; $i++ ) {
      my $site = join ":", $chr, $i;
      $hash{$site} = $e->Name;
    }
  }
  return \%hash;
}

# _processBedFile returns a matrix of the bedfile coordinates
sub _processBedFile {
  state $check = compile( Str, Str );
  my ( $self, $bedFile ) = $check->(@_);

  my @array;

  my @lines = path($bedFile)->lines( { chomp => 1 } );

  @lines = @lines[0 .. 400000 - 1];

  if(@lines > 400000) {
    $self->log('fatal', "Sorry, currently have limit of 400,000 ranges/primers");
  }

  # @lines = sort {$a->[0] <=> $}
  my $id = 0;
  for my $line (@lines) {
    my @fields = split /\t/, $line;
    my ( $chr, $start, $stop, $name ) = @fields;

    if ( defined $chr && defined $start && defined $stop ) {
      if ( !defined $name ) {
        $name = sprintf( "%s:%d-%d", $chr, $start, $stop );
      }

      $chr =~ s/\Achr//xmi;

      if ( $chr eq 'M' ) {
        $chr = 23;
      }
      elsif ( $chr eq 'X' ) {
        $chr = 24;
      }
      elsif ( $chr eq 'Y' ) {
        $chr = 25;
      }

      if(!looks_like_number($start) && looks_like_number($stop)) {
        $self->log('fatal', "chromStart and chromEnd must be numbers, found: $start, $stop"); 
      }
      
      $id++;

      if(length($name) > 31) {
        $name = substr($name, 0, 21) . "$id";
      }

      try {
        push @array, [$chr, $start, $stop, $name];
      }

      catch {
        $self->log( 'debug', "ignoring line: $line" );
      };

    }
    else {
      $self->log( 'fatal', sprintf( "Bedfile missing chr, start, or stop: %s", $line ) );
    }
  }
  return $self->_processBedObjs( \@array );
}

# _processBedObjs returns a matrix of the bedfile coordinates
sub _processBedObjs {
  state $check = compile( Str, ArrayRef );
  my ( $self, $bedObjAref ) = $check->(@_);

  my %sites;

  for my $b (@$bedObjAref) {
    my $chr = $b->[0];
    for ( my $i = $b->[1]; $i <= $b->[2]; $i++ ) {
      $sites{$chr}{$i} = $b->[3];
    }
  }

  undef @$bedObjAref;

  return $self->_bedSites( \%sites );
}

# _bedSites takes a hash of sites and creates a unique bedfile as an arrayref
sub _bedSites {
  state $check = compile( Str, HashRef );
  my ( $self, $sitesHref ) = $check->(@_);

  my ( %coveredSite, %coveredChr, @bed );
  my @chrs = ( 1 .. 26, 'M', 'X', 'Y' );

  for my $chr (@chrs) {
    if ( exists $sitesHref->{$chr} ) {
      $coveredChr{$chr} = 1;
      my @sites = sort { $a <=> $b } keys %{ $sitesHref->{$chr} };
      my @boundaries = $sites[0];
      for ( my $i = 1; $i < @sites; $i++ ) {
        if ( $sites[ $i - 1 ] != $sites[$i] - 1 ) {
          push @boundaries, $sites[ $i - 1 ], $sites[$i];
        }
      }
      push @boundaries, $sites[$#sites];

      # assign bed entries
      for ( my $i = 0; $i < @boundaries; $i += 2 ) {
        my $name  = $sitesHref->{$chr}{ $boundaries[$i] };

        # replicates
        # MPD::Bed::Raw->new(
        #   {
        #     Chr   => $chr,
        #     Start => $boundaries[$i],
        #     End   => $boundaries[ $i + 1 ],
        #     Name  => $name,
        #   }
        # );
        my $entry = [$chr, boundaries[$i], $boundaries[ $i + 1 ], $name];


        push @bed, $entry;

        # assign covered sites
        for ( my $pos = $boundaries[$i]; $pos < $boundaries[ $i + 1 ]; $pos++ ) {
          my $site = join ":", $chr, $pos;
          $coveredSite{$site} = $name;
        }
      }
    }
  }

  undef %$sitesHref;

  return ( \@bed, \%coveredChr, \%coveredSite );
}

# BUILDARGS takes either a hash reference or string, which is a primer file
sub BUILDARGS {
  my $class = shift;

  if ( scalar @_ == 1 ) {
    if ( !reftype( $_[0] ) ) {
      # assumption is that you passed a file be read and used to create
      # recall, reftype returns undef for string values

      my $file = $_[0];
      my ( $bedAref, $coveredChrHref, $coveredSitesHref ) = $class->_processBedFile($file);
      return $class->SUPER::BUILDARGS(
        {
          BedFile     => $file,
          Entries     => $bedAref,
          CoveredSite => $coveredSitesHref,
          CoveredChr  => $coveredChrHref,
        }
      );
    }
    elsif ( reftype( $_[0] ) eq "ARRAY" ) {
      my $bedObjsAref = $_[0];
      my ( $bedAref, $coveredChrHref, $coveredSitesHref ) =
        $class->_processBedObjs($bedObjsAref);
      return $class->SUPER::BUILDARGS(
        {
          Entries     => $bedAref,
          CoveredSite => $coveredSitesHref,
          CoveredChr  => $coveredChrHref,
        }
      );
    }
    elsif ( reftype( $_[0] ) eq "HASH" ) {
      return $class->SUPER::BUILDARGS( $_[0] );
    }
    else {
      return $class->log( 'fatal',
        "Construct MPD::Bed object with either a hashref or bed file" );
    }
  }
  else {
    return $class->log( 'fatal',
      "Construct MPD::Bed object with either a hashref or bed file" );
  }
}

__PACKAGE__->meta->make_immutable;

1;
