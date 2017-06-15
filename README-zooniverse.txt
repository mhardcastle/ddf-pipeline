Making image sets and uploading to Zooniverse

python ~/git/ddf-pipeline/fix_sourceid.py infile outfile 

(If source_id is broken)

In your image download directory:

python ~/git/ddf-pipeline/download_image_files.py file.fits

(Downloads all required files and makes the image list as file-list.txt.)

wc file-list.txt

(Find how many jobs you're going to be running)

qsub t 0-100 -v INFILE=file.fits lgz.qsub

python ~/git/ddf-pipeline/move_to_directory.py file-list.txt

Upload to Zooniverse...

export PANOPTES_PASSWORD=whatever
panoptes-subject-uploader ./manifest.csv --username mjh22 --project 2513 --workflow 1772
