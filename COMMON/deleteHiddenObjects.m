function deleteHiddenObjects(obj)
set(0,'showhiddenhandles','on')
delete(obj);
set(0,'showhiddenhandles','off')
