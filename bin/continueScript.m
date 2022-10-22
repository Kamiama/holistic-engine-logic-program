function[] = continueScript
% function that allows script to stop if user does not like current outputs

hold = 1;
next_step = input("\nContinue? (Y/N) >> ", 's');
while hold
    if upper(next_step) == "N"
        error("Program Ended")
        hold = 0;
    elseif upper(next_step) == "Y"
        hold = 0;
    else
        next_step = input("Please answer either Y or N >> ", 's');
    end
end