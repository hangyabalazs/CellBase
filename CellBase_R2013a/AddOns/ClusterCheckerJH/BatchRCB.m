DirName =uigetdir;

LIST=dir([DirName]); 

[s kk]=size(LIST);
for i=3:s
   name=char(LIST(i).name);

	if LIST(i).isdir==1;
	cd([DirName '\' name])
	RunClustBatch
	end

end