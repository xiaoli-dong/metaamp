#!/usr/bin/env perl
use Getopt::Long;
use strict;
use warnings;
use FindBin;

my ($taxontxt, $idir, $odir);
&GetOptions(
	    "i=s" =>\$taxontxt,
    "idir=s" =>\$idir,
    "odir=s" =>\$odir
	    
    );

($taxontxt) ||
    die "usage: $0 OPTIONS where options are: -i  <input file name: taxonmy summary file> -idir <input directory> -odir <output directory>\n";

my $bin = "$FindBin::RealBin";

my $jsonfile = "$taxontxt.json";
my $cmd = "$bin/taxons2jsontree.pl -i $idir/$taxontxt -d 3 > $odir/$jsonfile";
print STDERR "$cmd\n";
system($cmd);

my $taxontree = "$odir/$taxontxt.html";
print_taxonTree($taxontree);

my $taxonjs = "$odir/$taxontxt.js";
print_taxonjs($taxonjs, "$idir/$taxontxt", $odir);
chmod(0755, "$odir/$jsonfile") or die "could not chang $odir/$jsonfile permission, $!\n";
chmod(0755, $taxonjs) or die "could not chang $taxonjs permission, $!\n";



sub print_taxonTree{

    my ($taxontree) = @_;
    open (OUT, ">$taxontree") or die "Could not open $taxontree file to wrtie, $!\n";
    
    print OUT <<END;
    <html>
	<head>
        <meta http-equiv="Content-Type" content="text/html; charset=utf-8">
        <title>Taxonomy Tree</title>
	<link rel="stylesheet" type="text/css" href="http://ebg.ucalgary.ca/metaamp/css/ext-theme-classic-all.css">
	<script type="text/javascript" src="http://ebg.ucalgary.ca/metaamp/js/extjs/ext-all.js"></script>
	<script type="text/javascript" src="http://ebg.ucalgary.ca/metaamp/js/extjs/bootstrap.js"></script>
	<script type="text/javascript" src="$taxontxt.js"></script>
	
        <style type="text/css">
	body {font-size:100%;}
    h1 {font-size:2.0em;}
    h2 {font-size:1.875em;}
    h3 {font-size:1.15em;}
    p {font-size:0.875em;}
    
    #tree-example{
    
  background: #FFFFFF;
    //border: 2px solid #ccc;
	//margin-top: 50px;
    margin-left: 10px;
  padding: 10px 15px;
  color: #000;
    
}
#navigation {
 position: absolute;
 top: 0;
 left: 0;
 width: 10em;
}

#content {
margin-top: 10px;
margin-left: 10px;
padding: 10px 15px;
  color: #000;
}


</style>
    </head>
    <body>
    
    <div id="content">
    <h3>Taxonomic Tree</h3> 
    
    </div>
    <div id="loading-mask"></div>
    <div id="loading">
    <div class="loading-indicator">
    Loading...
    </div>
    </div>
    <div id="tree-example"></div>
    </body>
</html>
END
    
close(OUT);

}


sub print_taxonjs{
    
    my ($taxonjs, $taxontxt, $odir) = @_;
    open(INPUT, "$taxontxt") or die $!;
    #Taxonomy Total sample1 sample2
    my $head = <INPUT>;
    chomp $head;
    my @samples = ();
    if($head =~ /^#Taxonomy/){
	@samples = split(/\s+/, $head);
	
    }
    shift @samples;
    shift @samples;
    close(INPUT);
    
    open (OUT, ">$taxonjs") or die "Could not open $taxonjs file to wrtie, $!\n";
    
    my $entity = "var records = [";
    foreach my $sample (@samples){
	
	$entity .= "\n"
	    . "{\n"
	    ."data:{\n"
	    . "\"dataIndex\":\"$sample\",\n"
	    ." \"name\":\"$sample\",\n"
	    . "\"type\":\"float\"\n"
	    ."}\n"
	    ."},\n";
	
    }
    $entity =~ s/,$//;
    $entity .= "];";
    print OUT $entity;
    print OUT <<END;
    
    
    //I am building an application where the user can generate tables on the fly. The description of these tables are stored at the back-end.

//I needed to generate models from these descriptions and be able to use data grid objects to manipulate the data (crud operations). 

Ext.require([
    'Ext.data.*',
    'Ext.grid.*',
    'Ext.tree.*'
]);
var contextMenu = new Ext.menu.Menu({
	items: [{
		text: 'Retrieve Sequences',
		iconCls: 'edit'
		//handler: edit
	    },
	    {
		text: 'Hello',
		iconCls: 'edit'
		//handler: edit
	    },


	    ]
	
    });



// Skeleton store
var store_template = {
    proxy: {
        type: 'ajax',
        url: '$jsonfile'
	//url: '454AllContigs.fna.orf.nr.blat.json'
    },
    folderSort: true
};


myfunction = function(Mynode){ 
    Mynode.collapse(true);
    //alert("hh depth=" + depth);
    if(Mynode.hasChildNodes() == false || Mynode.getDepth() >= 2){ 
	return; 
	    }else if(Mynode.hasChildNodes()){ 
	Mynode.expand(false); 
	Mynode.eachChild(myfunction); 
    } 
} 
    
	//This Function recurses through Every Node till conditions are met
//Currently Depth is hard-coded to 4 but can be given a variable value to change as //per situation
    var gdepth = 2;
	myfunctionXiaoli = function(Mynode, depth){ 
	    Mynode.collapse(true);
	    if (typeof depth != 'undefined') {
		gdepth = depth;
	    }
	    
	    //alert("updaet gdepth=" + gdepth);
	    if(Mynode.hasChildNodes() == false || Mynode.getDepth() >= gdepth){ 
		 //alert("mynode=" + Mynode.getDepth());
		return; 
	    }else if(Mynode.hasChildNodes()){ 
		Mynode.expand(false); 
		Mynode.eachChild(myfunctionXiaoli); 
	    } 
	} 
	//this Function expands given node till depth given
	//As of now the depth is given 4 but can be changed as per requirement
	expandNodeTillDepth = function(treePanelId ,nodeId){
	    
	    var PNode = (Ext.getCmp(treePanelId)).getNodeById(nodeId); //this function is not working for me;
	    PNode.collapse(true);
	    PNode.expand(false);
	    PNode.eachChild(
			    myfunction
			    ); 
	}


var tree_template = {
    title: 'Taxonomy Tree Visualization',
    width: 1200,
    //height: 800,
    renderTo: 'tree-example',
    collapsible: true,
    useArrows: true,
    rootVisible: false,
    store: '_STORE_',
    multiSelect: true,
    //singleExpand: true,
    viewConfig: {
	stripeRows: true
    },
root: {
      
            expanded: true,
            loaded:true
        },	    
  
  
    dockedItems: '_DOCKEDITEMS_'
    
}
// Skeleton grid (_PLUGINS_ & _STORE_, are placeholders)
var grid_template = {
    columnWidth: 1,
    plugins: '_PLUGINS_',
    height: 300,
    store: '_STORE_'
}

// Skeleton window (_ITEMS_ is a placeholder)
var window_template = {
    title: 'Dynamic Model / Window',
    height: 400,
    width: 750,
    layout: 'fit',
    items: '_ITEMS_'
}
    var docked = '[{' + 
	    'xtype: \\'toolbar\\',' + 
	     'items: [{' + 
		'     text: \\'Expand All\\',' +
		 '   handler: function(){' +
			'tree.expandAll();' +
		    '}' +
		'}, {' +
		  '  text: \\'Collapse All\\',' +
		   ' handler: function(){' +
			'tree.collapseAll();' +
		    '}' +
		'},' +
		'\\'->\\',' +
		'\\'Expand the tree at level: \\',' +
		
    '{' +
	'xtype     : \\'combo\\',' +
	'width     : 100,' +
	'value: 4,' +
	'store     : [1,2,3,4,5,6,7],' +
	'listeners: {' +
	 '   select: function (combo, value) {' +
		'var depth = combo.getValue();' +
		
		'myfunctionXiaoli(tree.getRootNode(), depth);' +
	    '}' +
	'}' +
    '}	' +		    
	'	]' +
	'}]';
// Generate a model dynamically, provide fields
    function modelFactory(name, fields) {
	model =  {
	    extend: 'Ext.data.Model',
	    fields: fields
	};
	script = "Ext.define('"+name+"',"+Ext.encode(model)+");";
	//alert("model=" + script);
	eval(script);
	
    }


// Generate a dynamic store
function storeFactory(name,template,model){
    template.model = model;
    eval(name+" = Ext.create('Ext.data.Store',"+Ext.encode(template)+");");
}

// Generate a dynamic store
function treeStoreFactory(name,template,model){
    template.model = model;
    script = name+" = Ext.create('Ext.data.TreeStore',"+Ext.encode(template)+ ");";
    //alert("treeStore=" + script);
    eval(script);
    
}

// Generate a dynamic grid, .. store name is appended as a string because otherwise, Ext.encode
// will cause 'too much recursion' error (same for plugins)
function gridFactory(name,template,store,plugins){
    
    script =  name+" = Ext.create('Ext.grid.Panel', "+Ext.encode(template)+");";
    script = script.replace("\\"_STORE_\\"", store);
    script = script.replace("\\"_PLUGINS_\\"", plugins);
    eval(script);
}
function treeGridFactory(name,template,store){
    
    script =  name+" = Ext.create('Ext.tree.Panel',  "+Ext.encode(template)+ ");";
    script = script.replace("\\"_STORE_\\"", store);
    script = script.replace("\\"_DOCKEDITEMS_\\"", docked);
    //alert("treeGrid=" + script);
    eval(script);
    

}
// Generate a dynamic window, .. items are appended as a string to avoid Ext.encode error
function windowFactory(winName,winTemp,items){
    script += winName+" = Ext.create('Ext.window.Window',"+Ext.encode(winTemp)+").show();";
    script = script.replace("\\"_ITEMS_\\"", items);
    eval(script);
}

// Generate a model, a store a grid and a window dynamically from a record list!
function generateDynamicModel(records){
    
    fields = [
	      {name: 'term',type: 'string'}
	];

    columns = [{
	    xtype: 'treecolumn',
	    text: 'Taxonomic Lineage',
	    sortable: true,
	locked: true,
	width: 200,
	    dataIndex: 'term'
	    
	}];
    var count = 0;
    cm = [];
    for (var i = 0; i < records.length; i++) {
	 
        fields[i+1] =  {
            name: records[i].data.dataIndex,
            type: records[i].data.type
        };
	
        columns[i+1] = {
	    xtype: 'gridcolumn',
            text: records[i].data.name,
	    flex:1,
            sortable: true,
            dataIndex: records[i].data.dataIndex
            
        };
	cm[i+1] = records[i].data.dataIndex;
	count++;
    }
    fields[count+1] = {
	name: 'desc',
	type: "string"
    }
    tree_template.columns = columns;
   
    modelFactory('myModel',fields);
    treeStoreFactory('myStore',store_template,'myModel');
    treeGridFactory('tree',tree_template,'myStore');
    //windowFactory('myWindow',window_template,'[myGrid]');
    //alert("after make treeGrid");
    // Direct access to the store created above 
    myStore.load();
    tree.on('itemcontextmenu', function(view, record, item, index, event){
            //alert(record)
            //treePanelCurrentNode = record;
	    var nodesSelected = tree.getSelectionModel().selected.items;
	    if(nodesSelected.length > 0)
		{
		    var node = nodesSelected[0];
		    var desc = node.get('desc');
		    //alert('node.desc=' + desc);
		    uvar = "";
		    for (var i = 1; i < cm.length; i++) {
			
			tmp= cm[i] + '=' + node.get(cm[i]) + "&";
			uvar += tmp;
			
		    }
		    uvar +='desc=' + desc;
		    //alert(uvar);
		    var hh = record.get('term');
		    contextMenu.showAt(event.getXY());
		    contextMenu.on('click', function(contextMenu, item, term){
			    //var array = JSON.stringify([ 'foo', 'bar' ]);
			    var array = JSON.stringify(cm);
			    var url = 'http://example.com/?data=' + encodeURIComponent(array);
			    //alert('samples=' + url);
			    //alert('Item "' + item.id + '" at Menu: ' + contextMenu.id + ' was clicked.' + item.text + ' tree nodeid=' + Ext.encode(hh));
			    result = window.open("perl.cgi?" + uvar, "result","toolbar=1,width=900,height=600,top=0,left=0, status=1,location=1,menubar=1,scrollbars=1,resizable=1" );
			    result.focus();
			});
		}
            event.stopEvent();
	},this);

}

Ext.onReady(function(){
	//rowEditing = Ext.create('Ext.grid.plugin.RowEditing');
    
	setTimeout(function(){
		Ext.get('loading').remove();
		Ext.get('loading-mask').fadeOut({remove:true});
	    }, 350);
	
	//alert("hello");
	generateDynamicModel(records);
});

END
    
    close(OUT);

 
    
}
