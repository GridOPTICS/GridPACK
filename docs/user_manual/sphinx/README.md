To modify any parameters associated the read the docs documentation, start by
going to `readthedocs.org`. You should be able to do this by just typing
`readthedocs.org` into a browser.

Click on the `Log in` button and log in using your github credentials. Select
`Read the Docs Community`. Find GridPACK under your projects and click on the
link. Click on the `Settings` button. This will bring up a window with various
settings. You can edit any of these fields. When you are done, hit the `Save`
button at the bottom of the page.

#Old material from converting GridPACK documentation from LaTeX to restructured text

Convert GridPACK.tex to restructured text format by running

```
pandoc -f latex -t rst -s GridPACK.tex -o GridPACK.rst
```

in the `user_manual` directory.

To install sphinx: install the `sphinx-doc` and `sphinx-rtd-theme` on you
computer. On a Mac with homebrew, this can be done with
```
brew install sphinx-doc
```
```
brew install sphinx-rtd-theme
```
You may need to install additional packages using pip. If pip is not
available on your machine, you can install it using

```
python3 -m pip install
```

If you are using a virtual environment, you can do a pip
install in the environment. Otherwise, you may need to add a file
`.config/pip/pip.conf` in your home directory to enable pip to install
additional packages. Inside this file, add the lines

```
[global]
break-system-packages = true
```

Once sphinx is installed, cd into the docs/user_manual/sphinx directory and type

```
make html
```

Under the sphinx directory there will be a directory `_build/html`. Open the
`index.html` file in this directory use the `file://` syntax in a browser to
look at the web version of the documentation.
