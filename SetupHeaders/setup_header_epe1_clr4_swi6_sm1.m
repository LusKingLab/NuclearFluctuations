% Data names and colors.

nameNMBC = {'wt','epe1','clr4','swi6','swi6-sm1'};
nameMBC = {'wt_MBC','epe1_MBC','clr4_MBC','swi6_MBC','swi6-sm1_MBC'};
nameIntercleave = reshape([nameNMBC;nameMBC],1,2*length(nameNMBC));
nameAll = [nameNMBC,nameMBC];
nameAll2 = {'wild type','\itepe1\Delta','\itclr4\Delta','\itswi6\Delta',...
    '\itswi6-sm1','wild type MBC','\itepe1\Delta \rm MBC',...
    '\itclr4\Delta \rm MBC','\itswi6\Delta \rm MBC','\itswi6-sm1 \rm MBC'};
% colorblind barrier-free color pallet (Color Universal Design)
colorAll = [0,0,0;...        % wt
            102,45,145;...   % epe1
            247,148,29;...   % clr4
            0,165,115;...    % swi6
            218,111,171;...  % swi6-sm1
            0,0,0;...        % wt MBC
            102,45,145;...   % epe1 MBC
            247,148,29;...   % clr4 MBC
            0,165,115;...    % swi6 MBC
            218,111,171;]... % swi6-sm1 MBC
            /255;
numM = length(nameNMBC);
[ points2,faces2,psi,theta] = setup_img( );
p2um=0.16;
f2s=2.5;
[points,faces,edges,neighbors] = TriSphere(3);
