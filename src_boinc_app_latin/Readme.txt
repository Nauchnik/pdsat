CLIENT APP
[3.18] Added new version of encoding "minimum 72 cells"
	   MAX_NOF_RESTARTS_INC changed from 3010 to 4050
[3.03] 11.04.2013 [*] Minisat 2.2 version January 2013
[3.02] 02.11.2012 [*] fixed bug with unsaving results while checkpoint
[3.01] 11.09.2012 [*] MAX_NOF_RESTARTS_CLIENT changed from 8000 to 5000
[3.00] 02.09.2012 [*] diag10 istead of diag9. N = 10, rows_count = 3, K = 8. 26 literals in each of 20 tasks in WU
				  [+] POSITIVE_LITERALS_MIN_SIZE = 26. diag9_2 WUs will be interrupted, cause one's value is 23
[2.02] 04.08.2012 [*] ������������ ���� � ��� diag9_2 �������� ����� (���� ��������� ����� ���� ����������)
[2.01] 10.07.2012 [*] ���������� Ox, �������� SSE
[2.00] 23.06.2012 [+] minisat 2.2 instead of 14.1
[1.19] 04.06.2012 [*] to clear db in large diag tests: nof_learnts = solver_nclauses(s) / 3; - > nof_learnts   = solver_nclauses(s) / 8;
[1.16] 22.05.2012 [+] ������� ������� �� A5/1
				  [+] ��������� ��������
[1.15] 25.01.2011 [+] ������� ����� 50000 ���������� � ����� ��� 3 ��������� ��������� ����������
				  [*] ������ ��� ������ ������ (�������� valgrind)
{1.11] 14.11.2011 [*] DC_fractionDone - ��� � 256
[1.09] 13.10.2011 [+] �������� ����� �������� full_mask, part_mask � value
				  [*] DC_fractionDone - ��� � 1024
[1.07] 05.10.2011 [*] ����� �������� ������ ������� ���������� �� ������������ ������
			      [*] DC_fractionDone - ��� � 8192
[1.03] 23.09.2011 [+] �������� ����� ���. ������ � ���� c �����������
[1.02] 23.09.2011 [+] �������� checkpoint
				  [*] ������ ������ DC_fractionDone ���������� ��� � 4096 ������������ ����� �� range_val_count
[1.01] 17.09.2011 ��������� ������ DC_fractionDone
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