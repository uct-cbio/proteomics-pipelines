#!/bin/sh

if [ "$1" == "" ]
then
  echo "Please provide a config"
  echo
  exit
else
 config=`readlink -m $1`
  . $config
fi

project_name=$project_folder
project_dir=$project_dir
department=$department
user=$user
sub_dirs=$sub_dirs

echo "Start..."
for i in "${sub_dirs[@]}";
do
  echo "Rsyncing $project_dir/$i /mnt/researchdata/$department space mounted on dev-igisoro.cbio.uct.ac.za ..."
  cmd="rsync --size-only --append -Pavvvzh $project_dir/$i $user@dev-igisoro.cbio.uct.ac.za:/mnt/researchdata/$department/ProjectData/$project_name"
  echo $cmd
  eval $cmd
done

echo "Done."

