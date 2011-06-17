function saveDataFTP(filename)
    load ftplogin
    f = ftp('underwater.ece.ncsu.edu',ftpusername,ftppassword);
    cd(f,'/Volume_1/Personal/wccox')
    
    mput(f,filename);
end