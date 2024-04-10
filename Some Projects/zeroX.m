function [value,isterminal,direction] = zeroX(~,y)
    value = y(2);
    isterminal = 1;
    direction = 0;
end