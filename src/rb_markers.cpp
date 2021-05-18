#include <cstdio>
#include <getopt.h>
#include <string>
#include <chrono>
#include <thread>
#include <condition_variable>
#include <mutex>
#include <future>
extern "C" {
#include <zlib.h>
#ifndef AC_KSEQ_H
#include "kseq.h"
#include "kstring.h"
KSEQ_INIT(gzFile, gzread);
#endif
}
#include "rowbowt.hpp"
#include "rowbowt_io.hpp"
#include "thread_stream.hpp"

struct RbAlignArgs {
    std::string inpre = "";
    std::string fastq_fname = "";
    int ftab = 0;
    size_t wsize = 10;
    size_t max_range = 1000;
    size_t min_range = 0;
    size_t threads = 1;
};

void print_help() {
    fprintf(stderr, "rb_markers_only");
    fprintf(stderr, "Usage: rb_markers_only [options] <index_prefix> <input_fastq_name>\n");
    fprintf(stderr, "    --inexact                                     \n");
    fprintf(stderr, "    --wsize            <int>         window size for performing marker queries along read\n");
    fprintf(stderr, "    --max-range        <int>         range-size upper threshold for performing marker queries\n");
    fprintf(stderr, "    --min-range        <int>         range-size upper threshold for performing marker queries\n");
    fprintf(stderr, "    <input_prefix>                   index prefix\n");
    fprintf(stderr, "    <input_fastq>                    input fastq\n");
}

RbAlignArgs parse_args(int argc, char** argv) {
    int c;
    char* end;
    RbAlignArgs args;
    static struct option long_options[] {
        {"wsize", required_argument, 0, 'w'},
        {"max-range", required_argument, 0, 'r'},
        {"min-range", required_argument, 0, 'm'},
        {"threads", required_argument, 0, 't'}
    };
    int long_index = 0;
    while((c = getopt_long(argc, argv, "o:w:r:hft:", long_options, &long_index)) != -1) {
        switch (c) {
            case 't':
                args.threads = std::atol(optarg);
            case 'f':
                args.ftab = 1; break;
            case 'r':
                args.max_range = std::atol(optarg);
                break;
            case 'm':
                args.min_range = std::atol(optarg);
                break;
            case 'w':
                args.wsize = std::atol(optarg);
                break;
            case 'h':
                print_help();
                exit(0);
                break;
            default:
                print_help();
                exit(1);
                break;
        }
    }

    if (argc - optind < 2) {
        fprintf(stderr, "no argument provided\n");
        exit(1);
    }

    args.inpre = argv[optind++];
    args.fastq_fname = argv[optind++];
    return args;
}


uint8_t seq_ntoa_table[] = {
 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',
 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',
 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',
 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',
 'N', 'A', 'N', 'C', 'N', 'N', 'N', 'G', 'N', 'N', 'N', 'N', 'N', 'N', 'A', 'N',
 'N', 'N', 'N', 'N', 'T', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',
 'N', 'A', 'N', 'C', 'N', 'N', 'N', 'G', 'N', 'N', 'N', 'N', 'N', 'N', 'A', 'N',
 'N', 'N', 'N', 'N', 'T', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',
 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',
 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',
 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',
 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',
 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',
 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',
 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',
 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',
};


char comp_tab[] = {
    0,    1,    2,    3,    4,    5,    6,    7,    8,    9,    10,   11,   12,   13,   14,   15,
    16,   17,   18,   19,   20,   21,   22,   23,   24,   25,   26,   27,   28,   29,   30,   31,
    32,   33,   34,   35,   36,   37,   38,   39,   40,   41,   42,   43,   44,   45,   46,   47,
    48,   49,   50,   51,   52,   53,   54,   55,   56,   57,   58,   59,   60,   61,   62,   63,
    64,   'T',  'V',  'G',  'H',  'E',  'F',  'C',  'D',  'I',  'J',  'M',  'L',  'K',  'N',  'O',
    'P',  'Q',  'Y',  'S',  'A',  'A',  'B',  'W',  'X',  'R',  'Z',  91,   92,   93,   94,   95,
    64,   't',  'v',  'g',  'h',  'e',  'f',  'c',  'd',  'i',  'j',  'm',  'l',  'k',  'n',  'o',
    'p',  'q',  'y',  's',  'a',  'a',  'b',  'w',  'x',  'r',  'z',  123,  124,  125,  126,  127
};



struct KSeqString {
    void load(kseq_t* s) {
        seq.clear(); name.clear(); qual.clear();
        seq.assign(s->seq.s, s->seq.l);
        name.assign(s->name.s, s->name.l);
        qual.assign(s->qual.s, s->qual.l);
    }
    KSeqString() {}
    KSeqString(kseq_t* s) {
        load(s);
    }
    KSeqString& operator=(const KSeqString& rhs) {
        seq = rhs.seq;
        name = rhs.name;
        qual = rhs.qual;
        return *this;
    }
    std::string seq;
    std::string name;
    std::string qual;
};

/* code copied from seqtk
 * https://github.com/lh3/seqtk/blob/7c04ce7898ad5909bd309c6ba3cd9c3bd0651f0e/seqtk.c#L1464
 */
void revc_in_place(KSeqString& seq) {
    int c0, c1;
    for (size_t i = 0; i < seq.seq.size()>>1; ++i) { // reverse complement sequence
        c0 = comp_tab[(int)seq.seq[i]];
        c1 = comp_tab[(int)seq.seq[seq.seq.size() - 1 - i]];
        seq.seq[i] = c1;
        seq.seq[seq.seq.size() - 1 - i] = c0;
    }
    if (seq.seq.size() & 1) // complement the remaining base
        seq.seq[seq.seq.size()>>1] = comp_tab[(int)seq.seq[seq.seq.size()>>1]];
    if (seq.qual.size()) {
        for (size_t i = 0; i < seq.qual.size()>>1; ++i) // reverse quality
            c0 = seq.qual[i], seq.qual[i] = seq.qual[seq.qual.size() - 1 - i], seq.qual[seq.qual.size() - 1 - i] = c0;
    }
}

bool marker_cmp(MarkerT a, MarkerT b) {
    if (get_seq(a) == get_seq(b) && get_pos(a) == get_pos(b)) {
        return get_allele(a) < get_allele(b);
    } else if (get_seq(a) == get_seq(b)) {
        return get_pos(a) < get_pos(b);
    } else {
        return get_seq(a) < get_seq(b);
    }
}



class ThreadPool {
    public:
    ThreadPool(int n, const rbwt::RowBowt& r, const RbAlignArgs& a) 
        : nthreads(n)
        , rbwt(r)
        , args(a) 
    {
        for (int i = 0; i < n; ++i) {
            make_thread(i);
        }
    }
    ~ThreadPool() {
        stop = true;
        task_status.notify_all();
        for (auto& worker: workers) {
            worker.join();
        }
    }

    void add_task(kseq_t* kseq) {
        {
            std::unique_lock<std::mutex> lock(task_mutex);
            if (stop) {
                throw std::runtime_error("add_task stopped ThreadPool");
            }
            seq_queue.emplace(kseq);
        }
        task_status.notify_one();
    }

    private:
    void make_thread(int i) {
        auto worker = [this, i]() {
            std::vector<MarkerT> markers;
            KSeqString seq;
            std::ostringstream out_buf;
            int j = 0;
            bool rev;
            auto out_fn = [&](rbwt::RowBowt::range_t p, std::pair<size_t, size_t> q, std::vector<MarkerT> mbuf) {
                size_t qstart = rev ? seq.seq.size()-q.first-1 : q.first;
                out_buf << seq.name << " " << p.second-p.first+1 << " " << (rev ? "-" : "+") << " " << q.first << " " << q.second << " " << q.second-q.first+1;
                if (p.second - p.first + 1 >= this->args.min_range && mbuf.size()) {
                    for (auto m: mbuf) {
                        out_buf << " " << get_seq(m) << "/" <<  get_pos(m) << "/" << static_cast<int>(get_allele(m));
                        markers.push_back(m);
                    }
                }  else out_buf << " .";
                out_buf << std::endl;
            };
            while (true) {
                // let Pool know that a thread is free to take a task
                {
                    std::unique_lock<std::mutex> lock(this->task_mutex);
                    this->task_status.wait(lock, [this]() {return this->stop || !this->seq_queue.empty();});
                    if (this->stop && this->seq_queue.empty()) break;
                    // copy over data
                    seq = std::move(seq_queue.front());
                    seq_queue.pop();
                }
                markers.clear();
                // do work on data here
                for (auto& c: seq.seq) {
                    c = seq_ntoa_table[c];
                }
                rev = false;
                this->rbwt.get_markers_greedy_seeding(seq.seq, this->args.wsize, this->args.max_range, out_fn);
                revc_in_place(seq);
                rev = true;
                this->rbwt.get_markers_greedy_seeding(seq.seq, this->args.wsize, this->args.max_range, out_fn);
                j += 1;
                if (out_buf.tellp() > 4096) {
                    tout << out_buf.str();
                    out_buf.str(std::string());
                }
            }
            if (out_buf.tellp()) {
                tout << out_buf.str();
                out_buf.str(std::string());
            }
        };
        workers.push_back(std::thread(worker));
    }
    int nthreads=1;
    std::atomic_bool stop = false;
    std::mutex task_mutex;
    std::condition_variable task_status;
    std::condition_variable thread_status;
    std::vector<std::thread> workers;
    std::queue<KSeqString> seq_queue;
    const RbAlignArgs args;
    const rbwt::RowBowt& rbwt;
};



rbwt::RowBowt load_rbwt(const RbAlignArgs args) {
    auto flag = rbwt::LoadRbwtFlag::MA;
    if (args.ftab) {
        flag = flag | rbwt::LoadRbwtFlag::FT;
    }
    rbwt::RowBowt rbwt(rbwt::load_rowbowt(args.inpre, flag));
    return rbwt;
}



// TODO turn this into something that's compatible with threading
void rb_report(const rbwt::RowBowt& rbwt, const RbAlignArgs args, KSeqString seq, std::vector<MarkerT>& markers) {
    // get rid of Ns in seq
    for (size_t i = 0; i < seq.seq.size(); ++i) {
        seq.seq[i] = seq_ntoa_table[seq.seq[i]];
    }
    /*
    markers.clear();
    bool rev = false;
    rbwt.get_markers_greedy_seeding(seq.seq, args.wsize, args.max_range, fn);
    */
    auto out = tout;
    out << seq.seq << std::endl;
    revc_in_place(seq);
    bool rev = true;
    out << seq.seq << std::endl;
    /*
    rbwt.get_markers_greedy_seeding(seq.seq, args.wsize, args.max_range, fn);
    */
}

void rb_markers_all(RbAlignArgs args) {
    auto start = std::chrono::high_resolution_clock::now();
    rbwt::RowBowt rbwt(load_rbwt(args));
    auto stop = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> diff = stop - start;
    std::cerr << "loading rowbowt + markers took: " << diff.count() << " seconds\n";
    // initiate thread pool
    ThreadPool pool(args.threads, rbwt, args);
    // open file
    int err, nreads = 0, noccs = 0;
    gzFile fq_fp(gzopen(args.fastq_fname.data(), "r"));
    if (fq_fp == NULL) {
        fprintf(stderr, "invalid file\n");
        exit(1);
    }
    kseq_t* seq(kseq_init(fq_fp));
    // start assigning reads to threads
    int i = 0;
    while ((err = kseq_read(seq)) >= 0) {
        pool.add_task(seq);
        i += 1;
    }
    // error checking here
    switch(err) {
        case -2:
            fprintf(stderr, "ERROR: truncated quality string\n");
            exit(1); break;
        case -3:
            fprintf(stderr, "ERROR: error reading stream\n");
            exit(1); break;
        default:
            break;
    }
}

int main(int argc, char** argv) {
    rb_markers_all(parse_args(argc, argv));
}