
% plot velocity

function QC_01_12_CheckVelModel(dbpath, dbname, ID)

disp('By default reads r3huber catalogue from the database')

dbname_relocate = [dbpath 'relocate/relocate_use/' dbname '_r3huber'];

run('/opt/antelope/5.6/setup.m')

% open origin and associated phases database
db = dbopen(dbname_relocate, 'r');
db_assoc = dblookup(db, '', 'assoc', '', '');
db_arrival = dblookup(db, '', 'arrival', '', '');
db_origin = dblookup(db, '', 'origin', '', '');
db_AssocArrival = dbjoin(db_assoc, db_arrival);

% open origin table
db_origin = dbsubset(db_origin, ['orid==' num2str(ID)]);

EventTime = datenum(strtime(dbgetv(db_origin, 'time')));
AllOrids = dbgetv(db_origin, 'orid');
EventTime = EventTime(AllOrids == ID, :);

% get from the associated phase table
db_AssocArrival = dbsubset(db_AssocArrival, 'timedef=="d"');
db_AssocArrival = dbsubset(db_AssocArrival, ['orid==' num2str(ID)]);

AllAssocPha = dbgetv(db_AssocArrival, 'sta');
AllAssocTimeres = dbgetv(db_AssocArrival, 'timeres');
AllAssocIDs= dbgetv(db_AssocArrival, 'orid');
AllAssocDelta = dbgetv(db_AssocArrival, 'delta');
AllAssocTime = datenum(strtime(dbgetv(db_AssocArrival, 'time')));

% current assoc phases
CurrentAssocTime = AllAssocTime(AllAssocIDs == ID);
CurrentAssocDelta = AllAssocDelta(AllAssocIDs == ID);
CurrentAssocTimeres = AllAssocTimeres(AllAssocIDs == ID);

Time = (CurrentAssocTime - EventTime)*24*3600;
TheoTime = Time - CurrentAssocTimeres;
Distance = CurrentAssocDelta * 111.3;

mean(CurrentAssocTimeres)
median(CurrentAssocTimeres)

subplot(4,1,1:3)
plot(Distance, Time, '.')
hold on
plot(Distance, TheoTime, '-', 'Color', [.5 .5 .5])
xlabel('Distance [km]')
ylabel('Time since origin [s]')
xlim([0 200])
ylim([0 50])

grid on
subplot(4,1,4)
plot(Distance, CurrentAssocTimeres, '.')
xlabel('Distance [km]')
ylabel('Time Residual [s]')
hold on
grid on
xlim([0 200])

dbclose(db)








