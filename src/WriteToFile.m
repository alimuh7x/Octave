function WriteToFile(filename, varargin)
    % WRITE_VECTORS Write multiple equal-length vectors to a text file.
    % Usage:
    %   write_vectors("data.txt", x, y, z, ...)
    %
    % Each line contains elements from all vectors separated by tabs.

    n = numel(varargin);
    if n == 0
        error("Provide at least one vector.");
    end

    len = numel(varargin{1});
    for k = 2:n
        if numel(varargin{k}) ~= len
            error("All vectors must have the same length.");
        end
    end

    fid = fopen(filename, 'w');
    if fid == -1
        error("Cannot open file: %s", filename);
    end

    for i = 1:len
        for k = 1:n
            fprintf(fid, "%g\t", varargin{k}(i));
        end
        fprintf(fid, "\n");
    end

    fclose(fid);
    fprintf('Wrote %d vectors (%d rows) to %s\n', n, len, filename);
end

