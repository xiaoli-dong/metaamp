function defineAction(str){

     if(str=='otus'){
	 result = window.open("html/dummy.html", "result","toolbar=1,width=900,height=600,top=0,left=0, status=1,location=1,menubar=1,scrollbars=1,resizable=1" );
	 window.document.metaamp.target = "result";
	result.focus();
	window.document.metaamp.action="cgi-bin/run.cgi?action=otus";
    }
    if(str=='asv'){
	result = window.open("html/dummy.html", "result","toolbar=1,width=900,height=600,top=0,left=0, status=1,location=1,menubar=1,scrollbars=1,resizable=1" );
	window.document.metaamp.target = "result";
	result.focus();
	window.document.metaamp.action="cgi-bin/run.asv.cgi?action=asv";
    }

}

function displayHelp(str) {
    help = window.open(str, "help","toolbar=1,width=900,height=600,top=0,left=0, status=1,location=1,menubar=1,scrollbars=1,resizable=1");
    help.focus();


}
function displayDemo(str) {
    demo = window.open(str, "help","toolbar=1,width=900,height=600,top=0,left=0, status=1,location=1,menubar=1,scrollbars=1,resizable=1");
    demo.focus();


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

function help(divid, btnid, closeid){
    // Get the modal
    var modal = document.getElementById(divid);
    modal.style.display = "block";

    // Get the button that opens the modal
    var btn = document.getElementById(btnid);
    
    // Get the <span> element that closes the modal
    var span = document.getElementById(closeid);
    
    // When the user clicks the button, open the modal
    btn.onclick = function() {
	modal.style.display = "block";
    }
    
    // When the user clicks on <span> (x), close the modal
    span.onclick = function() {
	modal.style.display = "none";
    }
    
    // When the user clicks anywhere outside of the modal, close it
    window.onclick = function(event) {
	if (event.target == modal) {
	    modal.style.display = "none";
	}
    }
}
