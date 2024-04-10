function n = isnatural(t)
    if and (t-round(t) == 0,  t > 0)
        disp('n is natural');
    else 
        error('n is not natural');
    end
end

