To use the database hooks set the environment variable `DDF_PIPELINE_DATABASE` and make sure that
your home directory contains a file called `.surveys` with the database password in it.

Setting this variable has the following effects:

* download script will register a download of a file.
* unpack script will register unpacking
* make_mslist script will populate the database with details of the downloaded data
* pipeline will update status at the start and end of the run and on failure

You should also set `DDF_PIPELINE_CLUSTER` to the name of your cluster, e.g. 'Herts', 'Leiden'.

You can add fields to the database with status 'Preprocessed' with the `insert_field.py` script, e.g

```
insert_field.py 239680 1
```

(inserts field L239680 with priority 1 into the list). Automated or
manual queueing systems can then find the top-priority field that has
not yet been processed to download from SARA.
