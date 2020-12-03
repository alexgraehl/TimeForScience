
# Here's how to symlink these files into VSCode's location

# VSCode Settings directory
SDIR=$HOME/Library/Application\ Support/Code/User

mv -f $SDIR/settings.json $SDIR/settings.json.bak
mv -f $SDIR/keybindings.json $SDIR/keybindings.json.bak

ln -s $TIME_FOR_SCIENCE_DIR/Config/Alex_Williams/vscode/settings.json    $SDIR/
ln -s $TIME_FOR_SCIENCE_DIR/Config/Alex_Williams/vscode/keybindings.json $SDIR/

