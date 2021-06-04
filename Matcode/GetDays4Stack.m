function [timoffrot,tempoffs] = GetDays4Stack(fam)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function is to get the days and corresponding approximate zero
% crossings for stacking (all 7 stations), based on 40 sps, station order 
% is ['PGC 'SSIB ''SILB ' 'LZB  ' 'TWKB ' 'MGCB ' 'KLNB '];
% 
% when the plot shows that the dipole is later than 400, then plus the same
% amount to the tempoffs.
%
%
%
% Modified by Chao Song, chaosong@princeton.edu
% First created date:   2019/06/25
% Last modified date:   2019/06/25
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if isequal(fam,'002')  %DON'T FORGET additional family-specific delcarations around line 156
    
    % generate unique dates matrix that in that family
    timoffrot = Readbostock(fam);
    ind=[];
    timoffrot(ind, :)=[];
    
    % tempoffs are the zero-crossings at E component before rotation
    tempoffs=[1211 1297 1231 1221 1207 1186 1207]; %these are the zero crossings

elseif isequal(fam, '043')  % use 043 instead of 013, but with error
    
    % generate unique dates matrix that in that family
    timoffrot = Readbostock(fam);
    ind=[];
    timoffrot(ind, :)=[];
    
    tempoffs=[1290 1374 1349 1195 1210 1208 1254]; %these are the zero crossings
    
elseif isequal(fam, '141')  % use 141 instead of 025, but with error
        
    % generate unique dates matrix that in that family
    timoffrot = Readbostock(fam);
    ind=[];
    timoffrot(ind, :)=[];
    
    tempoffs=[1318 1396 1391 1182 1206 1212 1269]; %these are the zero crossings
    
elseif isequal(fam, '047')  % use 047 instead of 028
    
    % generate unique dates matrix that in that family
    timoffrot = Readbostock(fam);
    ind=[];
    timoffrot(ind, :)=[];
    
    tempoffs=[1240 1331 1274 1214 1213 1197 1227]; %these are the zero crossings
    
elseif isequal(fam, '010')  % use 010 instead of 056
    
    % generate unique dates matrix that in that family
    timoffrot = Readbostock(fam);
    ind=[10];
    timoffrot(ind, :)=[];
    
    tempoffs=[1288 1391 1296 1179 1207 1200 1203]; %these are the zero crossings
    
elseif isequal(fam, '144')  % use 144 instead of 084
    
    % generate unique dates matrix that in that family
    timoffrot = Readbostock(fam);
    ind=[];
    timoffrot(ind, :)=[];
    
    tempoffs=[1342 1411 1416 1191 1218 1226 1285]; %these are the zero crossings
    
    
elseif isequal(fam, '099')
    
    % generate unique dates matrix that in that family
    timoffrot = Readbostock(fam);
    ind=[8];
    timoffrot(ind, :)=[];
    
    tempoffs=[1315 1398 1378 1179 1206 1211 1268]; %these are the zero crossings
    
elseif isequal(fam, '068')  % use 068 instead of 115   
    
    % generate unique dates matrix that in that family
    timoffrot = Readbostock(fam);
    ind=[];
    timoffrot(ind, :)=[];
    
    tempoffs=[1281 1347 1352 1194 1202 1202 1242]; %these are the zero crossings
    
elseif isequal(fam, '125')
    
    % generate unique dates matrix that in that family
    timoffrot = Readbostock(fam);
    ind=[];
    timoffrot(ind, :)=[];
    
    tempoffs=[1262 1340 1316 1200 1205 1198 1235]; %these are the zero crossings for 002: PGC,SSIB,SILB.
    
elseif isequal(fam, '147')
    
    % generate unique dates matrix that in that family
    timoffrot = Readbostock(fam);
    ind=[];
    timoffrot(ind, :)=[];
    
    tempoffs=[1309 1379 1383 1191 1207 1212 1262]; %these are the zero crossings for 002: PGC,SSIB,SILB.
    
elseif isequal(fam, '017')    % use 017 instead of 149
    
    % generate unique dates matrix that in that family
    timoffrot = Readbostock(fam);
    ind=[2 7];
    timoffrot(ind, :)=[];
    
    tempoffs=[1311 1379 1389 1186 1204 1210 1262]; %these are the zero crossings for 002: PGC,SSIB,SILB.

    
% added 2020/02/12, for fam 243 with pgc trio 
elseif isequal(fam, '243')
    
    % generate unique dates matrix that in that family
    timoffrot= [2003 062; 
                2003 063;
                2004 196;
                2004 197;
                2004 198;
                2004 199;
                2005 254;
                2005 255;
                2005 256];
%     timoffrot = Readbostock(fam);
    
    tempoffs=[1241 1300 1270 1278 1262 1242 1238]; %these are the zero crossings

elseif isequal(fam, '001')    % use 017 instead of 149
    
    % generate unique dates matrix that in that family
    timoffrot = Readbostock(fam);
    ind=[];
    timoffrot(ind, :)=[];
    
    tempoffs=[1267 1393 1366 1171 1204 1208 1266]; %these are the zero crossings 
    
    
elseif isequal(fam, '019')  % use 047 instead of 028
    
    % generate unique dates matrix that in that family
    timoffrot = Readbostock(fam);
    ind=[];
    timoffrot(ind, :)=[];
    
    tempoffs=[1258 1344 1293 1194 1205 1196 1227]; %these are the zero crossings
    
elseif isequal(fam, '021')  % use 047 instead of 028
    
    % generate unique dates matrix that in that family
    timoffrot = Readbostock(fam);
    ind=[];
    timoffrot(ind, :)=[];
    
    tempoffs=[1260 1349 1300 1191 1204 1196 1231]; %these are the zero crossings 1253
    
elseif isequal(fam, '045')  % use 010 instead of 056
    
    % generate unique dates matrix that in that family
    timoffrot = Readbostock(fam);
    ind=[];
    timoffrot(ind, :)=[];
    
    tempoffs=[1300 1404 1355 1179 1218 1217 1289]; %these are the zero crossings    
    
elseif isequal(fam, '076')  % use 010 instead of 056
    
    % generate unique dates matrix that in that family
    timoffrot = Readbostock(fam);
    ind=[];
    timoffrot(ind, :)=[];
    
    tempoffs=[1300 1403 1365 1168 1203 1196 1258]; %these are the zero crossings  

elseif isequal(fam, '176')  % use 010 instead of 056
    
    % generate unique dates matrix that in that family
    timoffrot = Readbostock(fam);
    ind=[];
    timoffrot(ind, :)=[];
    
    tempoffs=[1295 1400 1353 1166 1200 1204 1254]; %these are the zero crossings
    
elseif isequal(fam, '015')  % use 141 instead of 025, but with error
        
    % generate unique dates matrix that in that family
    timoffrot = Readbostock(fam);
    ind=[];
    timoffrot(ind, :)=[];
    
    tempoffs=[1337 1412 1415 1170 1206 1218 1284]; %these are the zero crossings

elseif isequal(fam, '158')  % use 141 instead of 025, but with error
        
    % generate unique dates matrix that in that family
    timoffrot = Readbostock(fam);
    ind=[];
    timoffrot(ind, :)=[];
    
    tempoffs=[1331 1406 1410 1168 1207 1219 1290]; %these are the zero crossings
    
elseif isequal(fam, '234')  % use 141 instead of 025, but with error
        
    % generate unique dates matrix that in that family
    timoffrot = Readbostock(fam);
    ind=[];
    timoffrot(ind, :)=[];
    
    tempoffs=[1331 1406 1410 1168 1207 1219 1290]; %these are the zero crossings    
    tempoffs=tempoffs+148;

elseif isequal(fam, '231')  % use 141 instead of 025, but with error
        
    % generate unique dates matrix that in that family
    timoffrot = Readbostock(fam);
    ind=[];
    timoffrot(ind, :)=[];
    
    tempoffs=[1477 1552 1540 1297 1341 1356 1424]; %these are the zero crossings
    
elseif isequal(fam, '006')
    
    % generate unique dates matrix that in that family
    timoffrot = Readbostock(fam);
    ind=[];
    timoffrot(ind, :)=[];
    
    tempoffs=[1325 1399 1383 1174 1206 1213 1274]; %these are the zero crossings
    
else
    
    % generate unique dates matrix that in that family
    timoffrot = Readbostock(fam);
    timoffrot(ind, :)=[];
    
    tempoffs=[1315 1398 1378 1179 1206 1211 1268]; %these are the zero crossings

end

% 
% stas=['PGC  '
%     'SSIB '
%     'SILB '
%     'LZB  '
%     'TWKB '
%     'MGCB '
%     'KLNB '];















