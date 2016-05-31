### Process to setup backups
1. Get an account on dev-igosoro.cbio.uct.ac.za and asked to be added to the specific department/group. The /researchdata cifs mounts are on dev-igosoro.
2. Create ssh keypair on hexdata to be able to connect to dev-igosoro.
3. Each config file can have individual project parameter settings.
4. A cronjob needs to be added which calls the specific config e.g. 00 18 * * * /home/gerrit/projects/cbio-pipelines/16S/archive/archive_project.sh /home/gerrit/projects/cbio-pipelines/16S/archive/config.txt  2>&1 | mail -s "Cronjob report - Backup of immonology project test"  gerrit.botha@uct.ac.za

