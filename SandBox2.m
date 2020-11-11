
close all
clear all
clc


danny = Demo_Class("Danny", "Man");

john  = Demo_Class_handel("John", 'boy');
timmy = Demo_Class_handel("Timmy", 'girl');

addlistener(john,'NameChange',   @callback1);
addlistener(john,'MemberChange', @callback2);

addlistener(timmy,'NameChange',   @callback1);
addlistener(timmy,'MemberChange', @callback2);


tim  = timmy;
dave = danny;

dave.Change_Name("David");
john.Change_Name('Johana');
 tim.Change_Name("Tiger");

 
 
% Whoes who ??????
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


fprintf('\n\nTry Again......\n')
dave = dave.Change_Name("David");
fprintf('\nDave --> %s --> ', dave.name)
dave.Show_Member;
dave.Show_Number;


%% ======== Functions =====================================================
function callback1(src, eventData)

disp("Call Back -- Sombody changed their name")
disp(src.name)
disp(eventData.EventName)
src.Change_Number(rand);

end

function callback2(src, ~)

disp("Call Back -- Sombody changed their member!!")

src.Change_Member("Tranny")

end
