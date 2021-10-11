clear all
co = [];
for i = 1:100
    try
        General_case;
        co = [co,increase];
    end
end

result = mean(co);
dlmwrite('result_N',result,-append);