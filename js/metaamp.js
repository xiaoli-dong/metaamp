function defineAction(str){

    if(str=='rawstat'){
	//alert("rawstat");
	result = window.open("/metaamp/html/dummy.html", "result","toolbar=1,width=900,height=600,top=0,left=0, status=1,location=1,menubar=1,scrollbars=1,resizable=1" );
	window.document.metaamp.target = "result";
	result.focus();
	window.document.metaamp.action="/metaamp/cgi-bin/run.cgi?action=rawstat";
    }
    else if(str=='stripstat'){
	//alert("stripstat");
	result = window.open("/metaamp/html/dummy.html", "result","toolbar=1,width=900,height=600,top=0,left=0, status=1,location=1,menubar=1,scrollbars=1,resizable=1" );
	window.document.metaamp.target = "result";
	result.focus();
	window.document.metaamp.action="/metaamp/cgi-bin/run.cgi?action=stripstat";
    }

    else if(str=='otus'){
	//alert("otus");
	result = window.open("/metaamp/html/dummy.html", "result","toolbar=1,width=900,height=600,top=0,left=0, status=1,location=1,menubar=1,scrollbars=1,resizable=1" );
	window.document.metaamp.target = "result";
	result.focus();
	window.document.metaamp.action="/metaamp/cgi-bin/run.cgi?action=otus";
    }


}

function displayHelp(str) {
 help = window.open(str, "help","toolbar=1,width=900,height=600,top=0,left=0, status=1,location=1,menubar=1,scrollbars=1,resizable=1");
  help.focus();


}
function checkEmail(emailfield) {

    var myemail = document.getElementById(emailfield);
    var filter = /^([a-zA-Z0-9_\.\-])+\@(([a-zA-Z0-9\-])+\.)+([a-zA-Z0-9]{2,4})+$/;

    if (!filter.test(myemail.value)) {
    alert('Please provide a valid email address');
    myemail.focus;
    return false;
 }
    return true;
}
function checkAnalysisname(field){
   var name = document.getElementById(field);
   var filter = /^[A-Za-z0-9_]+/;
   if(!filter.test(name.value)){
     alert('Please provide a valide analysis name\nonly character, number and underscore are allowed');
     name.focus;
     return false;
   }
   return true;
}
function validateSubmit(){


    var x = document.forms["metaamp"]["analysisname"].value;
    if (x == null || x == "") {
        alert("Analysis Name must be filled out");
        return false;
    }
    checkAnalysisname("analysisname");
    var x = document.forms["metaamp"]["email"].value;
    if (x == null || x == "") {
        alert("Email address are required");
        return false;
    }
checkEmail("email");
    var x = document.forms["metaamp"]["seqArchiveFile"].value;
    if (x == null || x == "") {
        alert("Archived sequence file is required");
        return false;
    }
    var x = document.forms["metaamp"]["mappingFile"].value;
    if (x == null || x == "") {
        alert("sample2Sequences mapping file is required");
        return false;
    }


    return true;

}
