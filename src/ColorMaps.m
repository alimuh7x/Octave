classdef ColorMaps
    % ---------------------------------------------------------------------
    % ColorMaps: Enum-like container for standard colormap names
    %
    % Usage:
    %   colormap(ColorMaps.turbo);
    %   p.add_Image(x, y, data, "Title", ColorMaps.parula);
    % ---------------------------------------------------------------------

    properties (Constant)
        jet      = "jet";
        parula   = "parula";
        hot      = "hot";
        cool     = "cool";
        turbo    = "turbo";
        gray     = "gray";
        winter   = "winter";
        spring   = "spring";
        autumn   = "autumn";
        summer   = "summer";
        bone     = "bone";
        copper   = "copper";
        pink     = "pink";
        hsv      = "hsv";
        prism    = "prism";
    end
end

