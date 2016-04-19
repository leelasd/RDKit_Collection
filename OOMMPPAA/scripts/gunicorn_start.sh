#!/bin/bash
PROJECT_NAME=cocoa                                 # Name of the application
DJANGODIR=/home/vagrant/$PROJECT_NAME/src/WebApp        # Django project directory
SOCKFILE=/tmp/$PROJECT_NAME.sock  # we will communicte using this unix socket
USER=vagrant                                       # the user to run as
GROUP=vagrant                                    # the group to run as
NUM_WORKERS=5                                    # how many worker processes should Gunicorn spawn
DJANGO_SETTINGS_MODULE=WebApp.wsgisettings             # which settings file should Django use
DJANGO_WSGI_MODULE=WebApp.wsgi                     # WSGI module name
TIMEOUT=60
echo "Starting $PROJECT_NAME as `whoami`"


cd $DJANGODIR
# Activate the environment
export RDBASE=/srv/RDKit/rdkit
export LD_LIBRARY_PATH=$RDBASE/lib:$LD_LIBRARY_PATH
export DJANGO_SETTINGS_MODULE=$DJANGO_SETTINGS_MODULE
export PYTHONPATH=$DJANGODIR:$PYTHONPATH

# Create the run directory if it doesn't exist
RUNDIR=$(dirname $SOCKFILE)
test -d $RUNDIR || mkdir -p $RUNDIR

echo "testing"
# Start your Django Unicorn
# Programs meant to be run under supervisor should not daemonize themselves (do not use --daemon)
gunicorn ${DJANGO_WSGI_MODULE}:application \
  --name $PROJECT_NAME \
  --workers $NUM_WORKERS \
  --user=$USER --group=$GROUP \
  --log-level=debug \
  --bind=unix:$SOCKFILE\
  --timeout $TIMEOUT
