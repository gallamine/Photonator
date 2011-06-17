function loadStartFTP
    load ftplogin
    
    f = ftp('underwater.ece.ncsu.edu',ftpusername,ftppassword);
    cd(f,'/Volume_1/Personal/wccox')
    mget(f,'startupOptions.mat');
    
    
    close(f);
end