#! /usr/bin/perl
=head1 NAME
    
    obj2opengl - converts obj files to a VBO-friendly format
    
    =head1 SYNOPSIS
    
    obj2opengl [options] file
    
    use -help or -man for further information
    
    =head1 DESCRIPTION
    
    This script transforms a given OBJ file into a format suitable
    for OpenGL rendering. The output file format is easier to load
    than the OBJ format and can be encoded in binary to allow for
    faster load times.
    
    =head1 AUTHOR
    
    Marc Sunet (http://www.shellblade.net)
    
    Original work: Heiko Behrens (http://www.HeikoBehrens.net)
    
    =head1 VERSION
    
    14th July 2013
    
    =head1 COPYRIGHT
    
    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.
    
    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.
    
    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.

=head1 ACKNOWLEDGEMENTS

This script is based on the work of Margaret Geroch and Heiko Behrens.

=head1 REQUIRED ARGUMENTS

The first or the last argument has to be an OBJ file according
to this () specification.

=head1 OPTIONS

=over

=item B<-help>

Print a brief help message and exits.

=item B<-man>

Prints the extended manual page and exits.

=item B<-unit>

Fit the object in the unit cube, i.e. scale the object such the
longest dimension is 1 unit.

=item B<-scale <float>>

Sets the scale factor explicitly. Please be aware that negative numbers
are not handled correctly regarding the orientation of the normals.

=item B<-recenter>

Set the origin to the object's center.

=item B<-center <float> <float> <float>>

Set the object's center to the given point.

=item B<-o>, B<-outputFilename>

Name of the output file name. If omitted, the output file the same as the
input filename but with the extension .h

=item B<-nameOfObject>

Specifies the name of the generated variables. If omitted, same as
output filename without path and extension.

=item B<-noverbose>

Runs this script silently.

=item B<-bin>

Output in binary format.

=item B<-indexed>

Output indexed triangles.

=item B<-interleaved>

Output interleaved vertex attributes.

=item B<-matopt>

Optimise material output by encoding materials as vertex attributes.

=cut

use File::Basename;
use Getopt::Long;
use List::Util qw[min max];
use Pod::Usage;

# -----------------------------------------------------------------
# Main Program
# -----------------------------------------------------------------
handleArguments();

# derive center coords and scale factor if neither provided nor disabled
unless(defined($scalefac) && defined($xcen)) {
	calcSizeAndCenter();
}

if($verbose) {
	printInputAndOptions();
}

# TODO check integrity: Does every referenced vertex, normal and coord exist?
loadData();
normalizeNormals();

if($verbose) {
	printStatistics();
}

writeOutput();

# -----------------------------------------------------------------
# Sub Routines
# -----------------------------------------------------------------

sub handleArguments() {
    my $help = 0;
    my $man = 0;
    $unit = 0;
    $outputBinary = 0;
    $indexed = 0;
    $interleaved = 0;
    $matopt = 0;
    $verbose = 1;
    $errorInOptions = !GetOptions (
        "help" => \$help,
        "man"  => \$man,
        "unit" => \$unit,
        "scale=f" => \$scalefac,
        "center=f{3}" => \@center,
        "recenter" => \$recenter,
        "outputFilename=s" => \$outFilename,
        "nameOfObject=s" => \$object,
        "verbose!" => \$verbose,
        "bin" => \$outputBinary,
        "indexed" => \$indexed,
        "interleaved" => \$interleaved,
        "matopt" => \$matopt
        );
    
    if(!$unit) {
        $scalefac = 1;
    }
    
    if(defined(@center)) {
        $xcen = $center[0];
        $ycen = $center[1];
        $zcen = $center[2];
    }
    elsif (!$recenter) {
    	($xcen,$ycen,$zcen) = (0,0,0);
    }
    
    if($#ARGV == 0) {
        my ($file, $dir, $ext) = fileparse($ARGV[0], qr/\.[^.]*/);
        $inFilename = $dir . $file . $ext;
    } else {
        $errorInOptions = true;
    }
    
    # (optional) derive output filename from input filename
    unless($errorInOptions || defined($outFilename)) {
        my ($file, $dir, $ext) = fileparse($inFilename, qr/\.[^.]*/);
        $outFilename = $dir . $file;
        if ($outputBinary) {
            $outFilename = $outFilename . ".bin";
        }
        else {
            $outFilename = $outFilename . ".txt";
        }
    }
    
    # (optional) define object name from output filename
    unless($errorInOptions || defined($object)) {
        my ($file, $dir, $ext) = fileparse($outFilename, qr/\.[^.]*/);
        $object = $file;
    }
    
    ($inFilename ne $outFilename) or
        die ("Input filename must not be the same as output filename")
        unless($errorInOptions);
    
    if($errorInOptions || $man || $help) {
        pod2usage(-verbose => 2) if $man;
        pod2usage(-verbose => 1) if $help;
        pod2usage();
    }
    
    # check wheter file exists
    open ( INFILE, "<$inFilename" )
        || die "Can't find file '$inFilename' ...exiting \n";
    close(INFILE);
}

# Stores center of object in $xcen, $ycen, $zcen
# and calculates scaling factor $scalefac to limit max
#   side of object to 1.0 units
sub calcSizeAndCenter() {
    open ( INFILE, "<$inFilename" )
        || die "Can't find file $inFilename...exiting \n";
    
    $numVerts = 0;
    
    my ($xsum, $ysum, $zsum
       ,$xmin, $ymin, $zmin
       ,$xmax, $ymax, $zmax);
    
    while ( $line = <INFILE> ) {
        chop $line;
        if ($line =~ /v\s+.*/) {
            $numVerts++;
            @tokens = split(' ', $line);
            
            $xsum += $tokens[1];
            $ysum += $tokens[2];
            $zsum += $tokens[3];
            
            if ($numVerts == 1) {
                $xmin = $tokens[1];
                $xmax = $tokens[1];
                $ymin = $tokens[2];
                $ymax = $tokens[2];
                $zmin = $tokens[3];
                $zmax = $tokens[3];
			}
			else {
		        if ($tokens[1] < $xmin) {
		            $xmin = $tokens[1];
		        }
		        elsif ($tokens[1] > $xmax) {
		          $xmax = $tokens[1];
		        }
		        if ($tokens[2] < $ymin) {
		            $ymin = $tokens[2];
		       }
		        elsif ($tokens[2] > $ymax) {
		            $ymax = $tokens[2];
		        }
		        if ($tokens[3] < $zmin) {
		            $zmin = $tokens[3];
		        }
		        elsif ($tokens[3] > $zmax) {
		             $zmax = $tokens[3];
		         }
			}
        }
    }
    close INFILE;
    
    #  Calculate the center
    unless(defined($xcen)) {
        $xcen = $xsum / $numVerts;
        $ycen = $ysum / $numVerts;
        $zcen = $zsum / $numVerts;
    }
    
    #  Calculate the scale factor
    unless(defined($scalefac)) {
        my $xdiff = ($xmax - $xmin);
        my $ydiff = ($ymax - $ymin);
        my $zdiff = ($zmax - $zmin);

        if ( ( $xdiff >= $ydiff ) && ( $xdiff >= $zdiff ) ) {
            $scalefac = $xdiff;
        }
        elsif ( ( $ydiff >= $xdiff ) && ( $ydiff >= $zdiff ) ) {
            $scalefac = $ydiff;
        }
        else {
            $scalefac = $zdiff;
        }
        $scalefac = 1.0 / $scalefac;
    }
}

sub printInputAndOptions() {
    print "Input file     : $inFilename\n";
    print "Output file    : $outFilename\n";
    print "Object name    : $object\n";
    print "Center         : <$xcen, $ycen, $zcen>\n";
    print "Scale by       : $scalefac\n";
}

sub printStatistics() {
    my $numVerts = 3*$numFaces;
    my $ratio = $indexMapSize/$numVerts*100;
    print "----------------\n";
    print "Vertices       : $numVerts\n";
    print "Faces          : $numFaces\n";
    print "Texture Coords : $numTexture\n";
    print "Normals        : $numNormals\n";
    if ($indexed) {
        my $numIndices = scalar @indices;
        print "----------------\n";
        print "Indices                       : $numIndices\n";
        print "Output vertices (non indexed) : $numVerts\n";
        print "Output vertices (indexed)     : $indexMapSize\n";
        printf ("Size                          : %4.2f%\n", $ratio);
    }
    else {
        print "Vertices      : $numVerts\n";
    }
}

# reads vertices into $xcoords[], $ycoords[], $zcoords[]
#   where coordinates are moved and scaled according to
#   $xcen, $ycen, $zcen and $scalefac
# reads texture coords into $tx[], $ty[]
#   where y coordinate is mirrowed
# reads normals into $nx[], $ny[], $nz[]
#   but does not normalize, see normalizeNormals()
# reads faces and establishes lookup data where
#   va_idx[], vb_idx[], vc_idx[] for vertices
#   ta_idx[], tb_idx[], tc_idx[] for texture coords
#   na_idx[], nb_idx[], nc_idx[] for normals
#   store indizes for the former arrays respectively
#   also, $face_line[] store actual face string
# aabb left in $aabb
sub loadData {
    $numVerts = 0;
    $numFaces = 0;
    $numTexture = 0;
    $numNormals = 0;
    $numObjects = 0;
    $indexMapSize = 0;
    @indices = ();
    %vertexToIdx = ();
    
    open ( INFILE, "<$inFilename" )
        || die "Can't find file $inFilename...exiting \n";
    
    while ($line = <INFILE>) {
        chop $line;

        # vertices
        if ($line =~ /v\s+.*/) {
	        @tokens= split(' ', $line);
	        $x = ( $tokens[1] - $xcen ) * $scalefac;
	        $y = ( $tokens[2] - $ycen ) * $scalefac;
	        $z = ( $tokens[3] - $zcen ) * $scalefac;
	        $xcoords[$numVerts] = $x;
	        $ycoords[$numVerts] = $y;
	        $zcoords[$numVerts] = $z;
            
	        $numVerts++;
	        
	        if (defined($aabb)) {
	            $aabb->{min}->[0] = min($aabb->{min}->[0], $x);
	            $aabb->{min}->[1] = min($aabb->{min}->[1], $y);
	            $aabb->{min}->[2] = min($aabb->{min}->[2], $z);
                $aabb->{max}->[0] = max($aabb->{max}->[0], $x);
                $aabb->{max}->[1] = max($aabb->{max}->[1], $y);
                $aabb->{max}->[2] = max($aabb->{max}->[2], $z);
            }
            else {
                $aabb->{min} = [$x, $y, $z];
                $aabb->{max} = [$x, $y, $z];
            }
        }
        
        # texture coords
        elsif ($line =~ /vt\s+.*/) {
	    @tokens= split(' ', $line);
	    $x = $tokens[1];
	    $y = 1 - $tokens[2];
	    $tx[$numTexture] = $x;
	    $ty[$numTexture] = $y;
        
	    $numTexture++;
        }
        
        #normals
        elsif ($line =~ /vn\s+.*/) {
	    @tokens= split(' ', $line);
	    $x = $tokens[1];
	    $y = $tokens[2];
	    $z = $tokens[3];
	    $nx[$numNormals] = $x;
	    $ny[$numNormals] = $y;
	    $nz[$numNormals] = $z;
        
	    $numNormals++;
        }
        
        # faces
        elsif ($line =~ /f\s+([^ ]+)\s+([^ ]+)\s+([^ ]+)(\s+([^ ]+))?/) {
            @a = split('/', $1);
            @b = split('/', $2);
            @c = split('/', $3);
            $va_idx[$numFaces] = $a[0]-1;
            $ta_idx[$numFaces] = $a[1]-1;
            $na_idx[$numFaces] = $a[2]-1;
            
            $vb_idx[$numFaces] = $b[0]-1;
            $tb_idx[$numFaces] = $b[1]-1;
            $nb_idx[$numFaces] = $b[2]-1;
            
            $vc_idx[$numFaces] = $c[0]-1;
            $tc_idx[$numFaces] = $c[1]-1;
            $nc_idx[$numFaces] = $c[2]-1;
            
            $face_line[$numFaces] = $line;
            
            $numFaces++;
            
            # ractangle => second triangle
            if($5 != "") {
                @d = split('/', $5);
                $va_idx[$numFaces] = $a[0]-1;
                $ta_idx[$numFaces] = $a[1]-1;
                $na_idx[$numFaces] = $a[2]-1;
                
                $vb_idx[$numFaces] = $c[0]-1;
                $tb_idx[$numFaces] = $c[1]-1;
                $nb_idx[$numFaces] = $c[2]-1;
                
                $vc_idx[$numFaces] = $d[0]-1;
                $tc_idx[$numFaces] = $d[1]-1;
                $nc_idx[$numFaces] = $d[2]-1;
                
                $face_line[$numFaces] = $line;
                
                $numFaces++;
            }
            
            if ($indexed) {
                # Update vertex -> index map
                my $n = scalar keys %vertexToIdx;
                my $matIdx = $matIndexFromName{$material};
                my $vert1  = "$a[0]/$a[1]/$a[2]/$matIdx";
                my $vert2  = "$b[0]/$b[1]/$b[2]/$matIdx";
                my $vert3  = "$c[0]/$c[1]/$c[2]/$matIdx";
                my $vert4  = "";
                if ($5 != "") { $vert4 = "$d[0]/$d[1]/$d[2]/$matIdx"; }
                if (not exists $vertexToIdx{$vert1}) { $vertexToIdx{$vert1} = $n++; }
                if (not exists $vertexToIdx{$vert2}) { $vertexToIdx{$vert2} = $n++; }
                if (not exists $vertexToIdx{$vert3}) { $vertexToIdx{$vert3} = $n++; }
                if ($5 != "" and not exists $vertexToIdx{$vert4}) { $vertexToIdx{$vert4} = $n++; }
                
                # Push indices
                my $n = scalar @indices;
                $indices[$n++] = $vertexToIdx{$vert1};
                $indices[$n++] = $vertexToIdx{$vert2};
                $indices[$n++] = $vertexToIdx{$vert3};
                if ($5 != "")
                {
                    $indices[$n++] = $vertexToIdx{$vert1};
                    $indices[$n++] = $vertexToIdx{$vert3};
                    $indices[$n++] = $vertexToIdx{$vert4};
                }
            }
        }
        
        # groups
        elsif ($line =~ /usemtl\s*/) {
            @tokens = split(' ', $line);
            $material = $tokens[1];
            if ($indexed) {
                my $n = scalar @indices;
                $objects[$numObjects] = [$material, scalar @indices];
            }
            else {
                $objects[$numObjects] = [$material, $numFaces];
            }
            $numObjects++;
        }
        
        # material file
        elsif ($line =~ /mtllib/) {
            @tokens = split (' ', $line);
            $mtl = $tokens[1];
            loadMTL ($mtl);

            #print "Materials:\n";
            #for my $mat (@materials) {
            #	print "Material ".$mat->{"name"}."\n";
            #	print "Ka: ".$mat->{"Ka"}[0].", ".$mat->{"Ka"}[1].", ".$mat->{"Ka"}[2]."\n";
            #	print "Kd: ".$mat->{"Kd"}[0].", ".$mat->{"Kd"}[1].", ".$mat->{"Kd"}[2]."\n";
            #	print "Ks: ".$mat->{"Ks"}[0].", ".$mat->{"Ks"}[1].", ".$mat->{"Ks"}[2]."\n";
            #	print "Ns: ".$mat->{"Ns"}."\n";
            #}
        }
    }
    
    close INFILE;
    
    $indexMapSize = scalar keys %vertexToIdx;
    
    #print "Vertex -> Index:\n";
    #while (($key,$val) = each(%vertexToIdx)) { print "$key -> $val\n"; }
    #print "Indices:\n@indices\n";
}

sub normalizeNormals {
    for ( $j = 0; $j < $numNormals; ++$j) {
        $d = sqrt ( $nx[$j]*$nx[$j] + $ny[$j]*$ny[$j] + $nz[$j]*$nz[$j] );
        if ( $d == 0 ) {
	    $nx[$j] = 1;
	    $ny[$j] = 0;
	    $nz[$j] = 0;
        }
        else {
	    $nx[$j] = $nx[$j] / $d;
	    $ny[$j] = $ny[$j] / $d;
	    $nz[$j] = $nz[$j] / $d;
        }

    }
}

sub fixedIndex {
    local $idx = $_[0];
    local $num = $_[1];
    if($idx >= 0) {
        $idx;
    } else {
        $num + $idx + 1;
    }
}

sub loadMTL {
    my $i = -1;
    my $mtl = $_[0];
    $materials = ();
    
    print "Materials      : $mtl\n";
    
    open (MTL, "<$mtl")
        || die "Can't open file $mtl ... exiting\n";
    
    while (my $line = <MTL>) {
        chop $line;
        
        if ($line =~ /newmtl/) {
            @tokens = split (' ', $line);
            $mat = $tokens[1];
            $i++;
            $materials[$i] = {name => $mat};
        }
        elsif ($line =~ /^\s*Ns/) {
            @tokens = split (' ', $line);
            my $Ns = $tokens[1];
            $materials[$i]->{Ns} = $Ns;
        }
        elsif ($line =~ /^\s*Ka/) {
            @tokens = split (' ', $line);
            $materials[$i]->{Ka} = [$tokens[1], $tokens[2], $tokens[3]];
        }
        elsif ($line =~ /^\s*Kd/) {
            @tokens = split (' ', $line);
            $materials[$i]->{Kd} = [$tokens[1], $tokens[2], $tokens[3]];
        }
        elsif ($line =~ /^\s*Ks/) {
            @tokens = split (' ', $line);
            $materials[$i]->{Ks} = [$tokens[1], $tokens[2], $tokens[3]];
        }
    }
    
    close MTL;
    
    # Build mat index from name
    for (my $j = 0; $j < scalar @materials; $j++) {
        my $name = $materials[$j]{"name"};
        $matIndexFromName{$name} = $j;
    }
}

sub writeOutput {
    my $numMaterials = @materials;
    $matopt = $matopt && $numMaterials && $numObjects;
    my $writeMaterials = $numMaterials && !$matopt;
    my $writeObjects = $numObjects && !$matopt;
    
    open ( OUTFILE, ">:raw", $outFilename )
        || die "Can't create file $outFilename ... exiting\n";
    
    # write header
    print OUTFILE "verts";
    if ($numNormals)     { print OUTFILE " normals"; }
    if ($numTexture)     { print OUTFILE " texcoords"; }
    if ($outputBinary)   { print OUTFILE " binary"; }
    if ($indexed)        { print OUTFILE " indexed"; }
    if ($interleaved)    { print OUTFILE " interleaved"; }
    if ($writeMaterials) { print OUTFILE " materials"; }
    if ($writeObjects)   { print OUTFILE " objects"; }
    if ($matopt)         { print OUTFILE " matopt"; }
    print OUTFILE "\n";
    if ($indexed) {
        print OUTFILE "$indexMapSize\n";
        print OUTFILE scalar @indices."\n";
    }
    else {
        my $numVerts = 3*$numFaces;
        print OUTFILE "$numVerts\n"; # needed constant for glDrawArrays
    }
    
    # Write bounding box
    print OUTFILE "aabb $aabb->{min}->[0] $aabb->{min}->[1] $aabb->{min}->[2] ".
                       "$aabb->{max}->[0] $aabb->{max}->[1] $aabb->{max}->[2]\n";
    
    # materials
    if ($writeMaterials) {
        print OUTFILE "mats $numMaterials\n";
        for ($j = 0; $j < @materials; $j++) {
            my $name = $materials[$j]->{name};
            my $Ns   = $materials[$j]->{Ns};
            my $Ka   = $materials[$j]->{Ka};
            my $Kd   = $materials[$j]->{Kd};
            my $Ks   = $materials[$j]->{Ks};
            print OUTFILE "mat $name $Ns "
                ."$Ka->[0] $Ka->[1] $Ka->[2] "
                ."$Kd->[0] $Kd->[1] $Kd->[2] "
                ."$Ks->[0] $Ks->[1] $Ks->[2]\n";
        }
    }
    
    # compute object instances
    # objects are first sorted by material
    # vertex arrays or indices arrays are then reordered according to new object order
    # objects' offsets are then fixed
    # finally, objects sharing the same material are compressed into a single object
    if (!$matopt && $numObjects) {
        my $vertSize = 3;
        if ($numNormals) { $vertSize += 3; }
        if ($numTexture) { $vertSize += 2; }
        if ($matopt)     { $vertSize += 5; }
        for ($j = 0; $j < $numObjects; $j++) {
            my $matIndex = $matIndexFromName{$objects[$j][0]};
            if ($indexed) {
                $start = $objects[$j][1];
                if ($j < $numObjects-1) { $length = $objects[$j+1][1] - $start; }
                else { $length = scalar @indices - $start; }
            }
            else {
                if ($j < $numObjects-1) { $length = $objects[$j+1][1] - $objects[$j][1]; }
                else { $length = $numFaces - $objects[$j][1]; }
                $length = 3*$length;
                $start  = 3*$objects[$j][1];
                if ($interleaved)
                {
                    $start  *= $vertSize*4;
                    $length *= $vertSize*4;
                }
            }
            @objects[$j] = [$j, $matIndex, $start, $length];
        }
        
        # sort objects by material and group up
        @objects = sort { $a->[1] cmp $b->[1] } @objects;
        
        # re-order vertex arrays or indices
        if ($indexed) {
            my @indices_new = [];
            my $x = 0;
            for (my $j = 0; $j < scalar @objects; $j++) {
                my $start  = $objects[$j][2]; # Offset inside indices array
                my $length = $objects[$j][3]; # Length
                while ($length) {
                    $indices_new[$x++] = $indices[$start++];
                    $length--;
                }
            }
            @indices = @indices_new;
        }
        else {
            my @xcoords_new = [];
            my @ycoords_new = [];
            my @zcoords_new = [];
            my @nx_new = [];
            my @ny_new = [];
            my @nz_new = [];
            my @tx_new = [];
            my @ty_new = [];
            my $x = 0;
            for (my $j = 0; $j < scalar @objects; $j++) {
                my $start  = $objects[$j][2]; # Offset inside indices array
                my $length = $objects[$j][3]; # Length
                while ($length) {
                    $xcoords_new[$x] = $xcoords[$start];
                    $ycoords_new[$x] = $ycoords[$start];
                    $zcoords_new[$x] = $zcoords[$start];
                    $nx_new[$x] = $nx[$start];
                    $ny_new[$x] = $ny[$start];
                    $nz_new[$x] = $nz[$start];
                    $tx_new[$x] = $tx[$start];
                    $ty_new[$x] = $ty[$start];
                    $x++;
                    $start++;
                    $length--;
                }
            }
            @xcoords = @xcoords_new;
            @ycoords = @ycoords_new;
            @zcoords = @zcoords_new;
            @nx = @nx_new;
            @nz = @nz_new;
            @nz = @nz_new;
            @tx = @tx_new;
            @ty = @ty_new;
        }
        
        # recompute offsets
        my $n = scalar @objects;
        my $off = 0;
        for (my $j = 0; $j < $n; $j++) {
            $objects[$j][2] = $off;
            $off += $objects[$j][3];
        }
        
        # compress objects sharing the same material
        @objects_new = ();
        my $n = scalar @objects;
        my $mat = $objects[0][1];
        my $old_mat = $mat;
        my $start = 0;
        my $length = 0;
        for (my $j = 0; $j < $n; $j++) {
            #print "processing $objects[$j][0], $objects[$j][1], $objects[$j][2], $objects[$j][3]\n";
            $mat = $objects[$j][1];
            if ($mat != $old_mat) {
                #print "generating $old_mat, $start, $length\n";
                @objects_new[scalar @objects_new] = [0, $old_mat, $start, $length];
                $start += $length;
                $length = $objects[$j][3];
                $old_mat = $mat;
            }
            else {
                $length += $objects[$j][3];
            }
        }
        if ($mat == $old_mat) {
            #print "generating $old_mat, $start, $length\n";
            @objects_new[scalar @objects_new] = [0, $old_mat, $start, $length];
        }
        @objects = @objects_new;
        $numObjects = scalar @objects;
    }
    
    # object instances
    if ($writeObjects) {
        print OUTFILE "objs $numObjects\n";
        for (my $j = 0; $j < $numObjects; $j++) {
            print OUTFILE "obj $objects[$j][1] $objects[$j][2] $objects[$j][3]\n";
        }
    }
    
    # construct index to vertex map
    my @idxToVertex = ();
    if ($indexed) {
        while (my ($key,$val) = each(%vertexToIdx)) {
            $idxToVertex[$val] = $key;
        }
        for (my $i = 0; $i < scalar @idxToVertex; $i++) {
            #print "$i -> $idxToVertex[$i]\n";
        }
        #print "vertexToIdx keys: ".(scalar keys   %vertexToIdx)."\n";
        #print "vertexToIdx vals: ".(scalar values %vertexToIdx)."\n";
        #print "idxToVertex: ".(scalar @idxToVertex)."\n";
        #my $n = scalar keys %idxToVertex;
        #print "n = $n, indexMapSize = $indexMapSize\n";
        #print "Vertex -> Index:\n";
        #for (my $j = 0; $j < $n; $j++) { print "$j -> ".$idxToVertex{$j}."\n"; }
    }
    
    my $print = sub {
        my ($file, @data) = @_;
        if ($outputBinary) {
            my $n = scalar @data;
            print $file pack("f".$n, @data);
        }
        else { print $file join(" ", @data)."\n"; }
    };
    
    if ($interleaved) {
        if ($indexed) {
            for (my $j = 0; $j < $indexMapSize; $j++) {
                my $line = $idxToVertex[$j];
                #print "$j: line = $line\n";
                my @tokens = split ('/', $line);
                my $vi = fixedIndex ($tokens[0]-1, $indexMapSize);
                my $ni = -1;
                my $ti = -1;
                my $n = scalar @tokens;
                if ($numNormals) {
                    $ni = fixedIndex ($tokens[2]-1, $indexMapSize);
                }
                if ($numTexture) {
                    $ti = fixedIndex ($tokens[1]-1, $indexMapSize);
                }
                $print->(OUTFILE, $xcoords[$vi], $ycoords[$vi], $zcoords[$vi]); # verts
                if ($ni >= 0) { $print->(OUTFILE, $nx[$ni], $ny[$ni], $nz[$ni]); } # normals
                if ($ti >= 0) { $print->(OUTFILE, $tx[$ti], $ty[$ti]); } # tex coords
                if ($matopt)  { # mat attribs
                    my $mat = $materials[$tokens[$n-1]];
                    my $Kd  = $mat->{Kd};
                    my $Ks  = $mat->{Ks};
                    my $Ns  = $mat->{Ns};
                    my $K   = $materials[0]->{Kd};
                    #print "Kd: ".$K->[0].", ".$K->[1].", ".$K->[2]."\n";
                    $print->(OUTFILE, $Kd->[0], $Kd->[1], $Kd->[2], $Ks->[0], $Ns);
                }
            }
        }
        else {
            for(my $j = 0; $j < $numFaces; $j++) {
                $via = fixedIndex($va_idx[$j], $numVerts);
                $vib = fixedIndex($vb_idx[$j], $numVerts);
                $vic = fixedIndex($vc_idx[$j], $numVerts);
                if ($numNormals) {
                    $nia = fixedIndex($na_idx[$j], $numNormals);
                    $nib = fixedIndex($nb_idx[$j], $numNormals);
                    $nic = fixedIndex($nc_idx[$j], $numNormals);
                }
                if ($numTexture) {
                    $tia = fixedIndex($ta_idx[$j], $numTexture);
                    $tib = fixedIndex($tb_idx[$j], $numTexture);
                    $tic = fixedIndex($tc_idx[$j], $numTexture);
                }
                $print->(OUTFILE, $xcoords[$via], $ycoords[$via], $zcoords[$via]); # verts
                if ($ni >= 0) { $print->(OUTFILE, $nx[$nia], $ny[$nia], $nz[$nia]); } # normals
                if ($ti >= 0) { $print->(OUTFILE, $tx[$tia], $ty[$tia]); } # tex coords
                $print->(OUTFILE, $xcoords[$vib], $ycoords[$vib], $zcoords[$vib]); # verts
                if ($ni >= 0) { $print->(OUTFILE, $nx[$nib], $ny[$nib], $nz[$nib]); } # normals
                if ($ti >= 0) { $print->(OUTFILE, $tx[$tib], $ty[$tib]); } # tex coords
                $print->(OUTFILE, $xcoords[$vic], $ycoords[$vic], $zcoords[$vic]); # verts
                if ($ni >= 0) { $print->(OUTFILE, $nx[$nic], $ny[$nic], $nz[$nic]); } # normals
                if ($ti >= 0) { $print->(OUTFILE, $tx[$tic], $ty[$tic]); } # tex coords
            }
        }
    }
    # not interleaved
    else {
        if ($indexed) {
            for (my $j = 0; $j < $indexMapSize; $j++) {
                my $line = $idxToVertex[$j];
                @tokens = split ('/', $line);
                my $i = fixedIndex ($tokens[0]-1, $indexMapSize);
                $print->(OUTFILE, $xcoords[$i], $ycoords[$i], $zcoords[$i]);
                if ($numNormals) { $print->(OUTFILE, $nx[$i], $ny[$i], $nz[$i]); }
                if ($numTexture) { $print->(OUTFILE, $tx[$i], $ty[$i]); }
            }
        }
        else {
            for( $j = 0; $j < $numFaces; $j++) {
                $ia = fixedIndex($va_idx[$j], $numVerts);
                $ib = fixedIndex($vb_idx[$j], $numVerts);
                $ic = fixedIndex($vc_idx[$j], $numVerts);
                #print OUTFILE "\t\t// $face_line[$j]\n";
                $print->(OUTFILE, $xcoords[$ia], $ycoords[$ia], $zcoords[$ia]);
                $print->(OUTFILE, $xcoords[$ib], $ycoords[$ib], $zcoords[$ib]);
                $print->(OUTFILE, $xcoords[$ic], $ycoords[$ic], $zcoords[$ic]);
                if ($numNormals) {
                    $ia = fixedIndex($na_idx[$j], $numNormals);
                    $ib = fixedIndex($nb_idx[$j], $numNormals);
                    $ic = fixedIndex($nc_idx[$j], $numNormals);
                    $print->(OUTFILE, $nx[$ia], $ny[$ia], $nz[$ia]);
                    $print->(OUTFILE, $nx[$ib], $ny[$ib], $nz[$ib]);
                    $print->(OUTFILE, $nx[$ic], $ny[$ic], $nz[$ic]);
                }
                if ($numTexture) {
                    $ia = fixedIndex($ta_idx[$j], $numTexture);
                    $ib = fixedIndex($tb_idx[$j], $numTexture);
                    $ic = fixedIndex($tc_idx[$j], $numTexture);
                    $print->(OUTFILE, $nx[$ia], $ny[$ia]);
                    $print->(OUTFILE, $nx[$ib], $ny[$ib]);
                    $print->(OUTFILE, $nx[$ic], $ny[$ic]);
                }
            }
        }
    }
    
    # write indices
    if ($indexed) {
        for my $idx (@indices) {
            if ($outputBinary) { print OUTFILE pack("S<", $idx); }
            else { print OUTFILE "$idx\n"; }
        }
    }
    
    close OUTFILE;
}
