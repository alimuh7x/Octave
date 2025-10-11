function out = myprint(varargin)
    % myprintf("name1", val1, "name2", val2, ...)
    % Automatically formats:
    %  - small numbers (|val| < 1e4) as %.4f
    %  - large/small ones as %.3e
    %  - non-numeric as string
    out = "";
    for k = 1:2:length(varargin)
        name = varargin{k};
        val  = varargin{k+1};

        if isnumeric(val)
            if all(abs(val) < 1e4) && all(abs(val) > 1e-3)
                out = [out, sprintf('%s = %.4f  ', name, val)];
            else
                out = [out, sprintf('%s = %.3e  ', name, val)];
            end
        else
            out = [out, sprintf('%s = %s  ', name, num2str(val))];
        end
    end
    out = [out, '\n'];
    fprintf(out);
end
