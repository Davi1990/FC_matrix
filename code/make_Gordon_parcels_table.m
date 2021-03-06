function [T, cont, COMM_list, COMM_list_red] = make_Gordon_parcels_table(parcels)

%%
COMM_list = { % ordinamento mappe MC
    'Visual',...
    'RetrosplenialTemporal',...
    'SMhand',...
    'SMmouth',...
    'Auditory',...
    'CinguloOperc',...
    'VentralAttn',...
    'Salience',...
    'CinguloParietal',...
    'DorsalAttn',...
    'FrontoParietal',...
    'Default',...
    'None'...
};

COMM_list_red = { % ordinamento mappe MC
    'VIS',...
    'RSTP',...
    'SMH',...
    'SMM',...
    'AUD',...
    'CO',...
    'VAN',...
    'SAL',...
    'CP',...
    'DAN',...
    'FP',...
    'DMN',...
    'None'...
};


[num, txt, ~] = xlsread(parcels);
txtLabel      = txt(2:end,5);
txtHem        = txt(2:end,2);
real_pos      = [];
netID         = [];

for kk=1:13
    ROIid = cell2mat(COMM_list(kk));
    pos   = zeros(333,1);
    for AtlasNodePos=1:length(txtLabel)
        pos(AtlasNodePos)=strcmp(txtLabel{AtlasNodePos},ROIid);
    end
    pos      = find(pos);
    real_pos = [real_pos; pos];
    netID   = [netID; kk*ones(length(pos),1)];
    cont(kk) = length(pos);
end

for kk = 1: length(real_pos)
    HEMI_sorted{kk,1} = txtHem(real_pos(kk));
    NET_sorted{kk,1}  = txtLabel{real_pos(kk)};
end

T = table([1:333]',real_pos,NET_sorted,netID,HEMI_sorted,'VariableNames',{'ID' 'ID_sorted' 'NETName' 'NETID','HEM'});

