classdef Auxiliary
    methods (Static)
        function Interpreter(lang)
            arguments
                lang  {mustBeMember(lang,["latex","tex","none"])}
            end
            list_factory = fieldnames(get(groot,'factory'));
            index_interpreter = find(contains(list_factory,'Interpreter'));
            for i = 1:length(index_interpreter)
                default_name = strrep(list_factory{index_interpreter(i)},'factory','default');
                set(groot, default_name,lang);
            end
        end
    end
end