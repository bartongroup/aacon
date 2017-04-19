
$(document).ready(function(){
    $("#show_hidden1").click(function(){
        $("hidden_loader").show();
    });
    $("#show_hidden2").click(function(){
        $("hidden_loader").show();
    });
    $("#show_hidden3").click(function(){
        $("hidden_loader").show();
        $('html,body').scrollTop(0);
    });
    $("#show_hidden4").click(function(){
        $("hidden_loader").show();
        $('html,body').scrollTop(0);
    });
});