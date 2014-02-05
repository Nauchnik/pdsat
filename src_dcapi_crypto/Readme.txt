CLIENT APP
[4.00] 07.02.2014 [+] Added work with crypto problems
03.10.2013 - added using of known 2nd row for 1s square
           - in WU file fixed c ls_10_3_inc70
30.08.2013 - added generating of WUs with particular numbers (file with flag -w)
19.07.2013 - added generating of triple inc60
01.07.2013 - added writing problem_type (diag, inc0, inc70) into ls.problem_count
30.06.2013 - added generating of triple inc70
23.06.2013 - added new experiment: latin_triple
04.10.2013 - added limitaion in using of positive_literals; if values ended then stop creating new WUs
[3.03] 11.04.2013 [*] Minisat 2.2 version January 2013
01.03.2013 - interval of SQL launch was increased from 5 min to 1 hour
		   - if SQL fails with error, skip it and try again
[3.02] 02.11.2012 [*] fixed bug with unsaving results while checkpoint
[3.01] 11.09.2012 [*] MAX_NOF_RESTARTS_CLIENT changed from 8000 to 5000
[3.00] 02.09.2012 [*] diag10 istead of diag9. N = 10, rows_count = 3, K = 8. 26 literals in each of 20 tasks in WU
				  [+] POSITIVE_LITERALS_MIN_SIZE = 26. diag9_2 WUs will be interrupted, cause one's value is 23
[2.02] 04.08.2012 [*] Подгружаемый файл с КНФ diag9_2 уменьшен вдвое (были дубликаты почти всех дизъюнктов)
[2.01] 10.07.2012 [*] Компиляция Ox, включено SSE
[2.00] 23.06.2012 [+] minisat 2.2 instead of 14.1
[1.19] 04.06.2012 [*] to clear db in large diag tests: nof_learnts = solver_nclauses(s) / 3; - > nof_learnts   = solver_nclauses(s) / 8;
[1.16] 22.05.2012 [+] Пропуск заданий по A5/1
				  [+] Латинские квадраты
[1.15] 25.01.2011 [+] Останов после 50000 конфликтов и проба еще 3 вариантов начальной активности
				  [*] Убраны все утечки памяти (согласно valgrind)
{1.11] 14.11.2011 [*] DC_fractionDone - раз в 256
[1.09] 13.10.2011 [+] Добавлен вывод значений full_mask, part_mask и value
				  [*] DC_fractionDone - раз в 1024
[1.07] 05.10.2011 [*] Вывод значений только ядровых переменных из выполняющего набора
			      [*] DC_fractionDone - раз в 8192
[1.03] 23.09.2011 [+] Добавлен вывод вып. набора в файл c результатом
[1.02] 23.09.2011 [+] Добавлен checkpoint
				  [*] Вызовы вызовы DC_fractionDone происходят раз в 4096 срабатываний цикла по range_val_count
[1.01] 17.09.2011 Добавлены вызовы DC_fractionDone
[1.00] core_len -> mpi_base.core_len

MASTER APP
[1.01] 01/09/2012
[+] skip_values instead of skip_wus
[+] diag10 instead of diag9
[1.00] 11/01/2012
[+] parameter --results+to_skip= added for non creating new wus for some forst results
[*] adding new WU for every processed result, not only from current launch
[+] parameter --skip_wus=
[+] saving of WU_id
[*] FIRST_WU_CREATION changed to 65536
[+] procesing all CNFs in folder 
[*] for every CNF - own final_output file 