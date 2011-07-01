function saveDataFTP(filename,dir)
    load ftplogin
    f = ftp('underwater.ece.ncsu.edu',ftpusername,ftppassword);
    cd(f,'/Volume_1/Personal/wccox')
    cd(dir);
    mput(f,filename);
end