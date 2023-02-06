#
# Author: Miquel Anglada Girotto
# Contact: miquelangladagirotto [at] gmail [dot] com
#
# Script purpose
# --------------
# Download TCGA level3 data: gene expression - exon quantification corresponding
# to the samples in metadata
#
# Script arguments
# ----------------
#  - cancer_type
#  - download_dir
#  - data_category
#  - data_type
#  - workflow_type

require(optparse)
require(TCGAbiolinks)
require(tidyverse)
require(SummarizedExperiment)


##### FUNCTIONS #####
create_query = function(project_id, data_category, data_type, workflow_type){
    ## create query
    query = GDCquery(project       = project_id, 
                     data.category = data_category,
                     data.type     = data_type,
                     workflow.type = workflow_type)
    return(query)
}


download_data = function(project_id, download_dir, data_category, data_type, workflow_type){
    print(paste('Downloading...',project_id))
    
    query = create_query(project_id, data_category, data_type, workflow_type)
    GDCdownload(query, 
                method          = 'api', 
                directory       = download_dir,
                files.per.chunk = 10)
    print(paste('Finished downloading:',project_id))
}


parseargs = function(){
    
    option_list = list( 
        make_option("--cancer_type", type="character"),
        make_option("--data_category", type="character"),
        make_option("--data_type", type="character"),
        make_option("--workflow_type", type="character"),
        make_option("--download_dir", type="character")
    )

    args = parse_args(OptionParser(option_list=option_list))
    
    return(args)
}


main = function(){
    # process args
    args = parseargs()
    
    cancer_type = args[["cancer_type"]]
    project_id = paste('TCGA',cancer_type,sep='-')
    data_category = args[["data_category"]]
    data_type = args[["data_type"]]
    workflow_type = args[["workflow_type"]]
    download_dir = args[["download_dir"]]
    
    # run
    download_data(project_id, download_dir, data_category, data_type, workflow_type)
}

if(file.exists('MANIFEST.txt')) file.remove('MANIFEST.txt')


##### SCRIPT #####
if (sys.nframe() == 0L) {
    main()
    print("Done!")
}
