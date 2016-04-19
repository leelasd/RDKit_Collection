
PROJECT_NAME=$1

DB_NAME=$PROJECT_NAME
VIRTUALENV_NAME=$PROJECT_NAME

PROJECT_DIR=/home/vagrant/$PROJECT_NAME


cp -p $PROJECT_DIR/scripts/bashrc /home/vagrant/.bashrc


sudo apt-get install -y build-essential python-dev python-numpy git python-setuptools vim flex bison cmake sqlite3 libsqlite3-dev libboost-dev openbabel libboost-python-dev libboost-regex-dev python-matplotlib python-openbabel python-pip nginx libjpeg8 libjpeg62-dev libfreetype6 libfreetype6-dev unzip
sudo easy_install ipython Django==1.5.0 tornado
sudo pip install gunicorn django-jfu
sudo apt-get install -y python-psycopg2
sudo ln -s /usr/lib/x86_64-linux-gnu/libjpeg.so /usr/lib
sudo ln -s /usr/lib/x86_64-linux-gnu/libfreetype.so /usr/lib
sudo ln -s /usr/lib/x86_64-linux-gnu/libz.so /usr/lib
sudo pip install -U PIL

sudo mkdir /srv/RDKit
sudo chown -R  vagrant:vagrant /srv/RDKit
cd /srv/RDKit
git clone https://github.com/rdkit/rdkit.git
cd rdkit
git checkout dad1e1db5b0045e11972d931a8701c52ea90bb43

source ~/.bashrc

mkdir build
cd build
sudo cmake ..
sudo make
sudo make install

cd $RDBASE/External/INCHI-API
bash download-inchi.sh
cd $RDBASE/build
cmake .. -DRDK_BUILD_INCHI_SUPPORT=ON
make
make install


wget https://bitbucket.org/abradley/oommppaa/downloads/OOMMPPAA_db.db -O ~/$PROJECT_NAME/src/WebApp/data/OOMMPPAA_db.db

sudo apt-get install supervisor -y


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

sudo chown -R vagrant:vagrant /etc/nginx/sites*
echo 'upstream '$PROJECT_NAME'_server {
  # fail_timeout=0 means we always retry an upstream even if it failed
  # to return a good HTTP response (in case the Unicorn master nukes a
  # single worker for timing out).

  server unix:/tmp/'$PROJECT_NAME'.sock fail_timeout=0;
}

server {

    listen   80;
    server_name localhost;

    client_max_body_size 4G;

    access_log /tmp/nginx-access.log;
    error_log /tmp/nginx-error.log;

    location /static/ {
        alias   /home/vagrant/'$PROJECT_NAME'/src/WebApp/data/static/;
    }

    location /media/ {
        alias  /home/vagrant/'$PROJECT_NAME'/src/WebApp/data/data/media/;
    }

    location / {
        # an HTTP header important enough to have its own Wikipedia entry:
        #   http://en.wikipedia.org/wiki/X-Forwarded-For
        proxy_set_header X-Forwarded-For $proxy_add_x_forwarded_for;

        # enable this if and only if you use HTTPS, this helps Rack
        # set the proper protocol for doing redirects:
        # proxy_set_header X-Forwarded-Proto https;

        # pass the Host: header from the client right along so redirects
        # can be set properly within the Rack application
        proxy_set_header Host $http_host;

        # we dont want nginx trying to do something clever with
        # redirects, we set the Host: header above already.
        proxy_redirect off;

        # set "proxy_buffering off" *only* for Rainbows! when doing
        # Comet/long-poll stuff.  Its also safe to set if youre
        # using only serving fast clients with Unicorn + nginx.
        # Otherwise you _want_ nginx to buffer responses to slow
        # clients, really.
        # proxy_buffering off;

        if (!-f $request_filename) {
            proxy_pass http://'$PROJECT_NAME'_server;
            break;
        }
    }
}' >> /etc/nginx/sites-available/$PROJECT_NAME

sudo ln -s /etc/nginx/sites-available/$PROJECT_NAME /etc/nginx/sites-enabled/$PROJECT_NAME

sudo service nginx restart
sudo rm /etc/nginx/sites-enabled/default





