from subprocess import check_output
import re

def evaluate_params(delx, dely, eta, lt, theta, k,
	                kde_bw, kde_lags, kde_win):

    out = check_output(['Rscript', 'predict_with_params.r', str(delx), str(dely),
    	                str(eta), str(lt), str(theta), str(k),
    	                str(kde_bw), str(kde_lags), str(kde_win), 'fire'])
    
    #spearmint minimizes, so negate
    result = -float(re.sub(r'.*PEI: ([0-9.]*)', r'\1', out))
    
    print 'Result = %f' % result
    return result

def main(job_id, params):
    print params
    return evaluate_params(params['delx'], params['dely'], params['eta'], 
    	                   params['lt'], params['theta'], params['k'],
    	                   params['kde_bw'], params['kde_lags'], params['kde_win'])
