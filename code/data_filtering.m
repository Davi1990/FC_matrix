function new_data = data_filtering(data,hp,lp)

tr       = 2.6;
fs       = 1/tr; 
hp_new   = hp/(fs/2);
lp_new   = lp/(fs/2);
[b,a]    = butter(7,[hp_new lp_new],'bandpass');
new_data = filtfilt(b,a,data);