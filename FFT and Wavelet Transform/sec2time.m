function varargout = sec2time(sec)

HH = floor(sec/3600);
MM = mod(floor(sec/60), 60);
SSms = mod(sec,60);
if nargout <= 1
    varargout = cell(1,1);
    time = [HH MM SSms];
    varargout{1} = time;
else
    varargout = cell(3,1);
    varargout{1} = HH;
    varargout{2} = MM;
    varargout{3} = SSms;
end

end