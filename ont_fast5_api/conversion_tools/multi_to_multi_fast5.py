from __future__ import division

from argparse import ArgumentParser
from multiprocessing import Pool
from collections import deque
import logging
import h5py
import os

from ont_fast5_api import CURRENT_FAST5_VERSION, __version__
from ont_fast5_api.conversion_tools.conversion_utils import get_fast5_file_list, get_progress_bar
from ont_fast5_api.fast5_file import Fast5File
from ont_fast5_api.multi_fast5 import MultiFast5File

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)
exc_info = False


def batch_split_multi_to_adapter_multi(input_path, output_folder, threads, recursive, binning):

    pool = Pool(threads)
    file_list = get_fast5_file_list(input_path, recursive)
    pbar = get_progress_bar(len(file_list))

    def update(results):
        pbar.update(pbar.currval + 1)

    if not os.path.exists(output_folder):
        os.makedirs(output_folder)

    results_array = []
    for batch_num, filename in enumerate(file_list):
        results_array.append(pool.apply_async(multi_to_adapter_multi,
                                              args=(filename, output_folder, binning),
                                              callback=update))

    pool.close()
    pool.join()
    pbar.finish()

def multi_to_adapter_multi(input_file, output_folder, binning):
    results = deque([os.path.basename(input_file)])
    try:
        with MultiFast5File(input_file, 'r') as multi_f5:
            binned = {}
            for read_id in multi_f5.get_read_ids():
                if read_id in binning:
                    if binning[read_id] in binned:
                        binned[binning[read_id]].append(read_id)
                    else:
                        binned[binning[read_id]] = [read_id]
            for grp in binned:
                try:
                    output_file = os.path.join(output_folder, grp, os.path.basename(input_file))
                    if not os.path.exists(os.path.dirname(output_file)):
                        os.makedirs(os.path.dirname(output_file))
                    with MultiFast5File(output_file, 'a') as binned_multi_f5:
                        for read_id in binned[grp]:
                            if read_id not in binned_multi_f5.get_read_ids():
                                try:
                                    read = multi_f5.get_read(read_id)
                                    group_name = "read_" + read_id
                                    run_id = multi_f5.handle[group_name].attrs["run_id"]
                                    binned_read = binned_multi_f5.create_read(read_id, run_id)
                                    for group in read.handle:
                                        binned_read.handle.copy(read.handle[group], group)
                                except Exception as e:
                                    logger.error("{}\n\tFailed to add single read: '{}' to '{}'"
                                         "".format(e, read_id, output_file), exc_info=exc_info)
                except Exception as e:
                    logger.error("{}\n\tFailed to write to MultiRead file: {}".format(e, output_file), exc_info=exc_info)
            results.append(os.path.basename(input_file))
    except Exception as e:
        logger.error("{}\n\tFailed to copy files from: {}".format(e, input_file), exc_info=exc_info)
    finally:
        return results

def repackage_multi_to_multi(input_path, output_folder, recursive, filename_base, batch_size):
    file_list = get_fast5_file_list(input_path, recursive)
    pbar = get_progress_bar(len(file_list))
    def update():
        pbar.update(pbar.currval + 1)

    if not os.path.exists(output_folder):
        os.makedirs(output_folder)

    read_num = 0
    batch_num = 0
    output_file = os.path.join(output_folder, "{}_{}.fast5".format(filename_base, batch_num))
    output_multi_f5 = MultiFast5File(output_file, 'a')

    try:
        for input_file in file_list:
            with MultiFast5File(input_file, 'r') as multi_f5:
                for read_id in multi_f5.get_read_ids():
                    if read_num >= batch_size and batch_size > 0:
                        output_multi_f5.close()
                        batch_num += 1
                        output_file = os.path.join(output_folder, "{}_{}.fast5".format(filename_base, batch_num))
                        output_multi_f5 = MultiFast5File(output_file, 'a')
                        read_num = 0
                    try:
                        read = multi_f5.get_read(read_id)
                        group_name = "read_" + read_id
                        run_id = multi_f5.handle[group_name].attrs["run_id"]
                        output_read = output_multi_f5.create_read(read_id, run_id)
                        for group in read.handle:
                            output_read.handle.copy(read.handle[group], group)
                    except Exception as e:
                        logger.error("{}\n\tFailed to write to MultiRead file: {}"
                            "".format(e, output_file), exc_info=exc_info)
                    read_num += 1
            update()
    except:
        logger.error("{}\n\tFailed to join MultiRead files: {}"
            "".format(e, output_file), exc_info=exc_info)
    finally:
        output_multi_f5.close()
    print()
    return

def parse_csv(csv_fp):
    binning = {}
    with open(csv_fp, "r") as f:
        for line in f.readlines():
            read_id, grp = line.strip().split("\t")
            binning[read_id] = grp
    return binning

def main():
    parser = ArgumentParser("""Re-packages Multi-Fast5 files to Multi-Fast5 Files with the specified batch size.
                            Alternatively, each Multi-Fast5 File is split into multiple Multi-Fast5 files, each
                            containing reads of only one adapter group. The adapter group of each read is specified 
                            in a csv file supplied with argument '--binning_csv_file'. In the latter use case, 
                            the batch size argument is ignored.""")
    parser.add_argument('-i', '--input_path', required=True,
                        help="MultiRead fast5 file or path to directory of MultiRead files")
    parser.add_argument('-s', '--save_path', required=True,
                        help="Base directory to output folders of binned MultiRead fast5 files to")
    parser.add_argument('-b', '--binning_csv_file', required=False,
                        help="""CSV file (tab seperated, no header) containing read_ids in the 
                        first column and the corresponding binning group in the second column""")
    parser.add_argument('-f', '--filename_base', default='batch', required=False,
                        help="Root of output filename, default='batch' -> 'batch_0.fast5'")
    parser.add_argument('-n', '--batch_size', type=int, default=4000, required=False,
                        help="""Number of reads per multi-read file. If zero or a negative number is specified,
                         all reads are stored in a single file""")
    parser.add_argument('--recursive', action='store_true',
                        help="Search recursively through folders for MultiRead fast5 files")
    parser.add_argument('-t', '--threads', type=int, default=1, required=False,
                        help="Number of threads to use")
    parser.add_argument('-v', '--version', action='version', version=__version__)
    args = parser.parse_args()

    if args.binning_csv_file:
        binning = parse_csv(args.binning_csv_file)
        batch_split_multi_to_adapter_multi(args.input_path, args.save_path, args.threads, args.recursive, binning)
    else:
        repackage_multi_to_multi(args.input_path, args.save_path, args.recursive, args.filename_base, args.batch_size)


if __name__ == '__main__':
    main()