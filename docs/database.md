To use the database hooks set the environment variable `DDF_PIPELINE_DATABASE` and make sure that
your home directory contains a file called `.surveys` with the database password in it.

Setting this variable has the following effects:

* download script will register a download of a file.
* unpack script will populate the database with details of the downloaded data
* pipeline will update status at the start and end of the run

