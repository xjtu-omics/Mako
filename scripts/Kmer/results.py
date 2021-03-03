


def remove_repeat(hap_path):
    with open(hap_path) as fin:
        previous_line_str = ''
        previous_max_sim = -1
        previous_max_line = ''
        for line in fin.readlines():
            line_split = line.strip().split('\t')

            line_split[1] = str(int(line_split[1]) + 500)
            line_split[2] = str(int(line_split[2]) - 500)

            line_str = '-'.join(line_split[: 3])
            # print(line_str)

            if line_str != previous_line_str:
                print(previous_max_line)

                previous_line_str = line_str
                previous_max_sim = cal_overlap_ratio(int(line_split[1]), int(line_split[2]), int(line_split[3]), int(line_split[4]))
                previous_max_line = line.strip()
            else:
                if 'Valid' in line or 'Simple Event' in line:
                    sim = cal_overlap_ratio(int(line_split[1]), int(line_split[2]), int(line_split[3]), int(line_split[4]))
                    if sim > previous_max_sim:
                        previous_max_line = line.strip()
                else:
                    if previous_max_sim == -1:
                        previous_max_line = line.strip()

        print(previous_max_line)


def cal_overlap_ratio(start, end, valid_start, valid_end):
    size = end - start
    valid_size = valid_end - valid_start

    if size > valid_size:
        return round(valid_size / size, 2)
    else:
        return round(size / valid_size, 2)
if __name__ == '__main__':
    # hap1_path = '/mnt/e/Onedrive/stu/OneDrive - stu.xjtu.edu.cn/XJTU/Post/deepsv/results/valid/deepsv.s5.ssv/deepsv.s5.ssv.bed.vapor.na.revalid.hap1'
    # hap2_path = '/mnt/e/Onedrive/stu/OneDrive - stu.xjtu.edu.cn/XJTU/Post/deepsv/results/valid/deepsv.s5.ssv/deepsv.s5.ssv.bed.vapor.na.revalid.hap2'
    # remove_repeat(hap2_path)
    #
    # exit()
    #
    hap1_path = '/mnt/e/Onedrive/stu/OneDrive - stu.xjtu.edu.cn/XJTU/Post/mako/valid/NA19240/NA19240.mako_csvs.tsv.hap1'
    hap2_path = '/mnt/e/Onedrive/stu/OneDrive - stu.xjtu.edu.cn/XJTU/Post/mako/valid/NA19240/NA19240.mako_csvs.tsv.hap2'
    original_path = '/mnt/e/Onedrive/stu/OneDrive - stu.xjtu.edu.cn/XJTU/Post/mako/valid/NA19240/NA19240.mako_csvs.tsv'


    with open(original_path) as ori_file:

        valid_num = 0

        for line in ori_file:
            valid_flag1 = 0
            valid_sims = []

            line_split = line.strip().split('\t')

            chr = line_split[0]
            start = int(line_split[1])
            end = int(line_split[2])

            ex_start = start - 500
            ex_end = end + 500

            hap1_file = open(hap1_path)
            hap2_file = open(hap2_path)

            valid_info = []

            for line2 in hap1_file:
                line_split2 = line2.strip().split('\t')

                chr2 = line_split2[0]
                start2 = int(line_split2[1])
                end2 = int(line_split2[2])
                # print(chr, chr2, ex_start, start2, ex_end, end2)
                if chr == chr2 and ex_start == start2 and ex_end == end2:
                    valid_info.append(line_split2[3: ])
                    break
            for line2 in hap2_file:
                line_split2 = line2.strip().split('\t')

                chr2 = line_split2[0]
                start2 = int(line_split2[1])
                end2 = int(line_split2[2])

                if chr == chr2 and ex_start == start2 and ex_end == end2:
                    valid_info.append(line_split2[3:])
                    break


            for hap_valid_info in valid_info:

                if int(hap_valid_info[0]) > int(hap_valid_info[1]):
                    tmp = hap_valid_info[0]
                    hap_valid_info[0] = hap_valid_info[1]
                    hap_valid_info[1] = tmp

                if 'Valid' in hap_valid_info or 'Simple Event' in hap_valid_info:
                    size_sim = cal_overlap_ratio(start, end, int(hap_valid_info[0]), int(hap_valid_info[1]))
                    # print(size_sim)
                    hap_valid_info.append(str(size_sim))

                    valid_flag1 = 1
                    valid_sims.append(size_sim)
                else:
                    # pass
                    hap_valid_info.append('-1')
                    valid_sims.append(-1)


            if valid_flag1 == 1:
                valid_num += 1

            if valid_flag1 == 1:
                max_size_sim = valid_sims[0]
                max_size_sim_index = 0
                for i in range(1, len(valid_sims)):
                    sim = valid_sims[i]

                    if sim > max_size_sim:
                        max_size_sim = sim
                        max_size_sim_index = i
                # print(max_size_sim)

                print('{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\tValid'.format(chr, start, end, chr, valid_info[max_size_sim_index][0], valid_info[max_size_sim_index][1], max_size_sim), end='')
            else:
                pass
                if 'No reads in that region' in valid_info[0] and 'No reads in that region' in valid_info[1]:
                    print('{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\tNo reads'.format(chr, start, end, -1, -1, -1, -1), end='')
                else:
                    print('{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\tInValid'.format(chr, start, end, -1, -1, -1, -1), end='')

            # if valid_flag1 == 1:
            #     print(len(valid_sims))
            #     print('{0}\t{1}\t{2}\tValid\t{3}\t'.format(chr, start, end, max(valid_sims)), end='')
            # else:
            #     print('{0}\t{1}\t{2}\tInValid\t-1\t'.format(chr, start, end), end='')
            #
            # for hap_valid_info in valid_info:
            #     print('\t'.join(hap_valid_info), end='\t')


            print()
    # print(valid_num)


