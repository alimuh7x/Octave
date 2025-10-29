classdef StringBuilder < handle
    % StringBuilder: minimal, Octave-compatible class for building command text

    properties
        content = ""   % accumulated text
    end

    methods


        function append(obj, varargin)
            line = "";
            for i = 1:numel(varargin)
                value = varargin{i};
        
                if isnumeric(value)
                    part = num2str(value);
        
                elseif ischar(value) || isstring(value)
                    s = char(value);
        
                    % First argument (usually 'set', 'unset', 'plot') is command → no quotes
                    if i == 1 && (startsWith(s, "set") || startsWith(s, "unset") || startsWith(s, "plot"))
                        part = s;
                    else
                        % Always quote other string arguments
                        part = ["'", s, "'"];
                    end
        
                else
                    part = char(string(value));
                end
        
                line = [line, part, " "];
            end
        
            % Trim trailing space and add newline
            obj.content = [obj.content, strtrim(line), "\n"];
        end


        function clear(obj)
            obj.content = "";
        end

        function out = get(obj)
            out = char(obj.content);
        end

        function disp(obj)
            fprintf("%s", obj.content);
        end
    end
end
