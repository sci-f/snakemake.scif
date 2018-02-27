$(document).ready(function() {

    var github = "https://github.com/researchapps/living-docs"
    var divs = $("#h1,h2,h3,h4"); 
    $.each(divs, function(i,e){

        // Edit
        var did = $(e).attr('id');
        var start = "<div class='dropdown more'><span><i class='fa fa-ellipsis-h more' title='Edit'></i></span><div class='dropdown-content'>";

        // Permalink
        var link = "no-jekyll/#" + did;
        var button = "<p><a href='" + link + "' target='_blank'>Permalink</a></p>";
        start += button;

        // Edit
        var link = github + "/edit/master/docs/no-jekyll/#" + did;
        var button = "<p><a href='" + link + "' target='_blank'>Edit this page</a></p>";
        start += button;

        // Issues;
        var link = github + "/issues/new?labels=question&title=Question:&body=Question on: " + github + "/tree/master/docs/no-jekyll/README.md%23" + did;

        var button = "<p><a href='" + link + "' target='_blank'>Ask a Question</a></p>";
        start += button;
        start += "</div></div>";
        $(e).append(start)

    })
});
