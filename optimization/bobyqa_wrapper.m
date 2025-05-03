function bobyqa_est = bobyqa_wrapper(URA, Un, ref, wavelength, ...
                                     pos0, x_bounds, y_bounds, ...
                                     tol, max_iter)
%BOBYQA_WRAPPER  Interface MATLAB→Python (retorna sempre double 1×2).

    % (1) coloca dir. atual no sys.path (se necessário)
    if count(py.sys.path, pwd)==0
        insert(py.sys.path,int32(0),pwd);
    end

    % (2) importa e recarrega módulo
    driver = py.importlib.import_module('bobyqa_driver');
    driver = py.importlib.reload(driver);

    % (3) extrai limites reais se vetores completos forem passados
    xmin = min(x_bounds(:));  xmax = max(x_bounds(:));
    ymin = min(y_bounds(:));  ymax = max(y_bounds(:));

    % (4) listas simples para Python
    pos0_py = py.list({pos0(1),       pos0(2)});
    low_py  = py.list({xmin,          ymin});
    up_py   = py.list({xmax,          ymax});

    % (5) chamada (retorna ndarray float64)
    x_py = driver.bobyqa_optimize( ...
               py.numpy.array(URA), ...
               py.numpy.array(Un),  ...
               py.numpy.array(ref), ...
               wavelength, ...
               pos0_py, low_py, up_py, ...
               tol, int32(max_iter) );

    % (6) conversão imediata → elimina warnings de serialização
    bobyqa_est = double(x_py);   % garante vetor 1×2
end
