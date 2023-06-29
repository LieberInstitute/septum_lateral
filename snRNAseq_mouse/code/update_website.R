library("here")
library("withr")

## To clone septum_lateral_website we used:
# git clone git@github.com:LieberInstitute/septum_lateral.git --branch gh-pages --single-branch septum_lateral_website

with_dir(here(), {
    rmarkdown::render("README.Rmd", "html_document")
    system("mv README.html ../septum_lateral_website/index.html")
})

with_dir(
    gsub("septum_lateral", "septum_lateral_website", here()), {
    system("git commit -am -'Updated website with snRNAseq_mouse/code/update_website.R'; git push origin gh-pages")
})
