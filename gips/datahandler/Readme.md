## Installation

Currently the datahandler xmlrpc and logging daemon are started from an upstart script:

usage: `sudo service geokitd start`

`/etc/init/geokitd.conf`

```
# start geokitd 
# shawn patti

start on runlevel [2345]

script

export VIRTUAL_ENV="/usr/local/dh_venv"
export PATH="$VIRTUAL_ENV/bin:$PATH"

alias pydoc="python -m pydoc"
echo $$ > /var/run/geokitd.pid

exec /usr/local/dh_venv/bin/geokitd.py >> /var/log/geokitd.log 2>&1

end script

pre-start script
  echo "[`date`] Starting Geokitd" >> /var/log/geokitd.log
end script

pre-stop script
  rm /var/run/geokitd.pid
  echo "[`date`] Stopping Geokitd" >> /var/log/geokitd.log
end script
```
