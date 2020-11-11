
close all
clear all
clc


% --- Make Classes ----
danny = Demo_Class("Danny", "Man");

john  = Demo_Class_handel("John", 'boy');
timmy = Demo_Class_handel("Timmy", 'girl');

% --- Add listeners ----
addlistener(john,'NameChange',           @callback1);
addlistener(john,'MemberChange',         @callback2);
addlistener(john,'ObjectBeingDestroyed', @callback3);


addlistener(timmy,'NameChange',           @callback1);
addlistener(timmy,'MemberChange',         @callback2);
addlistener(timmy,'ObjectBeingDestroyed', @callback3);

% --- Copy classes ----
tim  = timmy;
dave = danny;

% --- Change Names ----
dave.Change_Name("David");
john.Change_Name('Johana');
 tim.Change_Name("Tiger");

 
% --- Randomize Time :( ---
tim = tim.rando;
 
 
% --- Whoes who ?????? ---
fprintf('\n\n\tWhos Who??????\n\n')

fprintf('\nDave --> %s --> ', dave.name)
dave.Show_Member;
dave.Show_Number;

fprintf('\nDanny --> %s --> ', danny.name)
danny.Show_Member;
danny.Show_Number;

fprintf('\nJohn --> %s --> ', john.name)
john.Show_Member;
john.Show_Number;

fprintf('\nTimmy --> %s --> ', timmy.name)
timmy.Show_Member;
timmy.Show_Number;

fprintf('\nTim --> %s --> ', tim.name)
tim.Show_Member;
tim.Show_Number;


% --- Lets try dave again ---
fprintf('\n\nTry Again......\n')
dave = dave.Change_Name("David");
fprintf('\nDave --> %s --> ', dave.name)
dave.Show_Member;
dave.Show_Number;

fprintf('\nDanny --> %s --> ', danny.name)
danny.Show_Member;
danny.Show_Number;

clear john

try
    delete(timmy)
catch E
    disp(E)
    clear timmy
end


%% Do it in a loop
clc

% --- Make people ----
count = 4;
for id = { "John", "Timmy", "Suzzy", "Andy";
           'boy',  'girl',  'girl',  'boy'}
    
    people(1:count) = Demo_Class_handel(id{1}, id{2});
    count = count - 1;
end

% --- Show who they are ----
for person = people
    fprintf('\nPerson: %s --> ', person.name)
    person.Show_Member;
    person.Show_Number;
end


% --- change them ----
newID = [ "Jannet", "Tess", "Hank", "Alese"; 
             " Boy",   "Boy",  "Girl",  "Girl"];

ind = 1;
for person = people
    
    addlistener(person,'NameChange',   @callback1);
    addlistener(person,'MemberChange', @callback2);
    
    person.Change_Identity(newID(1,ind), newID(2,ind));
    ind = ind + 1;
    
end


% --- Show who they've become
fprintf('\n\n\n')

for person = people
    fprintf('\nPerson: %s --> ', person.name)
    person.Show_Member;
    person.Show_Number;
end



%% ======== Functions =====================================================
function callback1(src, eventData)
fprintf('\n\n')
disp("Call Back -- Sombody changed their name")
disp(src.name)
disp(eventData.EventName)
src.Change_Number(69);
fprintf('\n')
end


function callback2(src, ~)

fprintf('\n\n')
disp("Call Back -- Sombody changed their member!!")

src.Change_Member("Tranny");
fprintf('\n')
end



function callback3(src, eventData)

fprintf('\n\n')
disp("Call Back -- Sombodys being destroyed!!!")
disp(src.name)
disp(eventData.EventName)
fprintf('\n')

end

