# NImGen

## Requirements:

* pandas
* abagen
* selenium (only for webgestalt)

## Installing

```
git clone https://github.com/juaml/nimgen.git
cd nimgen
python setup.py develop
```


## Using webgestalt

To use webgestalt, the respective selenium driver must be installed.

Firefox drivers: https://sites.google.com/a/chromium.org/chromedriver/downloads
Firefox: https://github.com/mozilla/geckodriver/releases

Just download the respective file, unzip and place in the users PATH (`/bin` or `/usr/bin`)

Note for Macos users: 
Open a Terminal, copy the driver to `/usr/local/bin` and then execute. An error message will appear. Press close or cancel. Open system preferences, security and press the button to allow for the execution. Go to the terminal, execute again. Another warning will appear. Allow.