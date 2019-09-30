#!/usr/bin/env perl
use Getopt::Long;
use strict;

my ($input);
&GetOptions("i=s"       => \$input
    );

($input) ||
    die "usage: $0 OPTIONS

where options are: -i  <newick format tree>\n";


#((29:0.0573123,27:0.0573123):0.0363852,28:0.0936975):0.406303;
#rooted tree is not working well with jsphylosvg, convert it to unrooted tree
open(INPUT, "$input") or die "Could not open $input to read, $!\n";
my $data = "";
while (<INPUT>){
    chomp;
    $data .= $_;
    
}
close(INPUT);
$data =~ s/:?\d+?\.\d+?;$/;/;
#print $data, "\n";


print  <<END;
<html>
    <head>
    <script type="text/javascript" src="../js/jsphylosvg-1.55/raphael-min.js" ></script> 
    <script type="text/javascript" src="../js/jsphylosvg-1.55/jsphylosvg-min.js"></script> 
    <script type="text/javascript" src="../js/jquery-2.1.3.min.js"></script>
    <script type="text/javascript">
    \$(document).ready(function(){
	//	$.get("tree.newick", function(data) {
	    var data = '$data';
	    
	    
	    var dataObject = {
	      newick: data,
	      fileSource: false
	    };		
	    phylocanvas = new Smits.PhyloCanvas(
		dataObject,
		'svgCanvas', 
		800, 800
		//'circular'
		);
	    
	
	    
	    //	});
	
		       });
</script>
    
    </head>
    
    
    <body>
  <div id="svgCanvas"> </div>
</body>
</html>		


END
