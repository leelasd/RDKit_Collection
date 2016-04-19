
PROJECT_NAME=$1

DB_NAME=$PROJECT_NAME
VIRTUALENV_NAME=$PROJECT_NAME

PROJECT_DIR=/home/vagrant/$PROJECT_NAME


cp -p $PROJECT_DIR/scripts/bashrc /home/vagrant/.bashrc



sudo chown -R  vagrant:vagrant /etc/supervisor/conf.d/
sudo chmod +x /home/vagrant/cocoa/scripts/gunicorn_start.sh
echo '[program:'$PROJECT_NAME'] 
command = bash '$PROJECT_DIR'/scripts/gunicorn_start.sh '$PROJECT_NAME'                 ; Command to start app
user = vagrant                                                          ; User to run as
stdout_logfile = /tmp/gunicorn_supervisor.log   ; Where to write log messages
redirect_stderr = true                                                ; Save stderr in the same log
environment=LANG=en_US.UTF-8,LC_ALL=en_US.UTF-8                       ; Set UTF-8 as default encoding' >> /etc/supervisor/conf.d/$PROJECT_NAME.conf

sudo supervisorctl reread
sudo supervisorctl update









