run_container () {
    # Load the Docker image from the tar.gz file
    docker pull ghcr.io/rgcgithub/regenie/regenie:v3.2.7.gz

    docker run \
      -e HOME=${WD} \
      -v ${WD}/:$HOME/ \
      ghcr.io/rgcgithub/regenie/regenie:v3.2.7.gz $cmd
}