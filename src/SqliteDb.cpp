/*
    SqliteDb.cpp
 *
 *  Created on: Sep 26, 2012
 *      Author: dima
 */
#include "SqliteDb.h"

#define BLOB_TYPE unsigned short int
#define NA 18446744073709551615

using namespace std;

int SqliteDb::intermediate_commit(){
//    xct->commit();
    db->execute("COMMIT");
//    delete [] xct;
//    xct = new sqlite3pp::transaction(*db, true, IMMEDIATE);
    db->execute(IMMEDIATE ? "BEGIN IMMEDIATE" : "BEGIN");
    return 0;
}

int SqliteDb::commit(){
     try { db->execute("COMMIT");}
     catch (exception& ex) { cerr << ex.what() << endl; return 1; }
     return 0;
}

void SqliteDb::new_transaction(){
    clog << "immediate: " << IMMEDIATE << endl;
    db->execute(IMMEDIATE ? "BEGIN IMMEDIATE" : "BEGIN");
    // xct = new sqlite3pp::transaction(*db, true, IMMEDIATE);
}
///////////////////////////////////////////////////////////////////////////////////
void SqliteDb::read_register_table(){
// one needs to check whether the table `register` exists whatsoever here:

//    const char * qrycheck = "SELECT  count(*) FROM sqlite_master WHERE type='table' AND name='register';";
//    sqlite3pp::query qry(*db, qrycheck);
    try {
        sqlite3pp::query qry(*db, "SELECT name, notes  FROM register");

        for (int i = 0; i < qry.column_count(); ++i) {
          cout << qry.column_name(i) << "\t";
        }
    } catch (int n) {
        clog << "`register` table was not found" << endl;
    }
}
void SqliteDb::parseSQLite() {
    return; }
///////////////////////////////////////////////////////////////////////////////////
void SqliteDb::init_sql_table( std::string & sample_label_ ){
    init_sql_table( sample_label_.c_str()  );
}

void SqliteDb::init_sql_table(const char * sample_label_  ){
    counter = 0;
    sample_label = std::string( sample_label_ );
    table_name = sample_label + "__vcf";

    char sql[256];
    sprintf(sql, "CREATE TABLE %s ("  \
                       "contig TEXT NOT NULL, " \
                       "pos INT  NOT NULL , " \
                       "ref TEXT, " \
                       "alt TEXT, " \
                       "id TEXT ," \
                       "qual FLOAT , "\
                       "filter TEXT );"\
                      , table_name.c_str() );
    exec_sql_log(sql);
}

void SqliteDb::init_register_table( ){
    const char * sql = "CREATE TABLE IF NOT EXISTS register ( "
                 "name TEXT PRIMARY KEY, "
            "vcf_file TEXT, "
            "bam_file TEXT, "
            "mapper TEXT, "
            "nm INT, "
            "mq INT, "
            "nt_q INT, "
            "wt BOOL, "
            "notes TEXT ); ";
   exec_sql_log(sql);
}

void SqliteDb::place_register_record_begin( const std::string & snpfile_ ,  const std::string & read_filename){
    snpfile = snpfile_;
    std::string note =  "coverage analysis started : " + snpfile;
    sqlite3pp::command cmd(*db, "INSERT OR REPLACE INTO register (name, vcf_file,bam_file, notes) VALUES (:name, :vcf_file, :bam_file, :notes) ");
    cmd.bind(":name", sample_label.c_str());
    cmd.bind(":vcf_file", snpfile.c_str());
    cmd.bind(":bam_file", read_filename.c_str());
    cmd.bind(":notes", note.c_str() );

// sqlite3pp::command cmd( *db, "INSERT OR REPLACE INTO register (name, vcf_file, notes) VALUES (:name, :vcf_file, :notes) ");
//    cmd.bind(":vcf_file", snpfile.c_str() );
    try {
        cmd.execute();
    } catch (exception& ex) {
        cerr << ex.what() << endl;
    }

}
void SqliteDb::place_register_record( ){
    std::string note =  "coverage_done : " + snpfile;
  {
     sqlite3pp::command cmd( *db, "INSERT OR IGNORE INTO register (name, notes) VALUES (:name, :notes) ");
    cmd.bind(":name", sample_label.c_str());
    // cmd.bind(":vcf_file", snpfile.c_str() );
    cmd.bind(":notes", note.c_str() );

    try {
        cmd.execute();
    } catch (exception& ex) {
        cerr << ex.what() << endl;
   }
  }
  {
    sqlite3pp::command cmd( *db, "UPDATE register SET notes = (:notes) WHERE name = (:name) " );
    cmd.bind(":name", sample_label.c_str());
    cmd.bind(":notes", note.c_str() );

    try {
        cmd.execute();
    } catch (exception& ex) {
        cerr << ex.what() << endl;
    }
  }
}

void SqliteDb::init_contig_table(){
    if (verbose){ clog << "contigs table:" << sample_label << endl;};

    char sql[256];
    sprintf(sql, "CREATE TABLE %s__contigs ("  \
                      "chr_name TEXT PRIMARY KEY     NOT NULL , " \
                       "data_type TEXT, " \
                       "chr_id INT, " \
                       "max_pos INT, " \
                       "records INT);", sample_label.c_str() );
    exec_sql_log(sql);
}

void SqliteDb::exec_sql_log(char const * sql){
    try {
        int exitcode = db->execute(sql); // parent class method
        if (verbose){   clog << "[sqlite3:] " << sql << " [returned:]  " << exitcode << endl;}
     } catch (exception& ex) {
         cerr << ex.what() << endl;
    }
}
/*
void SqliteDb::place_contig_table_record( size_t & chr_ref, string & chr_name ){

    char sql[256];
    sprintf(sql, "INSERT OR REPLACE INTO %s__contigs " \
        "(chr_name, data_type, chr_id) " \
        //"(chr_name, data_type, chr_id, records ) " 
        "VALUES " \
        "(:chr_name, :data_type, :chr_id);", sample_label.c_str() );
 
    sqlite3pp::command cmd( *db, sql);
    cmd.bind(":chr_name", chr_name.c_str() );
    cmd.bind(":data_type", "coverage");
    cmd.bind(":chr_id", (int) chr_ref + 1 );
//    cmd.bind(":records", n_snp );
    try {
        cmd.execute();
    } catch (exception& ex) {
        cerr << ex.what() << endl;
    }
}
*/

void SqliteDb::composite_index(){
    char sql[256];
    sprintf( sql, "CREATE UNIQUE INDEX idx1 ON %s (contig, pos);", table_name.c_str() );
    exec_sql_log(sql);
}

void SqliteDb::print_cov_db( const vcf::Variant & var ){
                        // clog << "printing ..." << endl;
//                    cout << "REACHED 1!" << endl;
/*
    size_t nbytes = sizeof(BLOB_TYPE);

    size_t buf_len = cov.range*2 + 1 ; // class int 
    size_t buf_tot_len = 4 * buf_len; // 4 is the number of subarrays, independent of BLOB_TYPE 

    BLOB_TYPE aln_ctr_str[buf_tot_len] ;
    cov.aln_centres->sprint_subarrays(aln_ctr_str );
*/
    counter ++;
    char sql[1024];
    sprintf (sql, "INSERT OR REPLACE INTO %s "  \
         "(contig, pos, id, ref, alt, qual, filter) " \
         "VALUES (:contig, :pos, :id, :ref, :alt, :qual, :filter) ", \
          table_name.c_str() );

    std::cout<< "contig " << var.sequenceName << "\tpos: " <<  var.position << endl;

    sqlite3pp::command cmd( *db, sql);
    cmd.bind(":contig", var.sequenceName.c_str() );
    cmd.bind(":pos", (const long long int) var.position );
    cmd.bind(":id", var.id.c_str() );
    cmd.bind(":ref", var.ref.c_str() );
    // report the list of alternate alleles.
/*    std::stringstream ss;
    std::ostream buff_ea(ss.rdbuf());
    var.printAlt(buff_ea);
    std::string buff_i;
    buff_i << ss.str();
    cmd.bind(":alt", buff_i.c_str() );
    */
    cmd.bind(":qual", var.quality);
    cmd.bind(":filter", var.filter.c_str() );

//   cmd.bind(":aln_ctr", (void const*) aln_ctr_str, buf_tot_len * nbytes);

    try {
        cmd.execute();
    } catch (exception& ex) {
        cerr << ex.what() << endl;
    }

    if (! counter % commit_every){
          intermediate_commit();  
        }
}

std::string & stream2string(std::ostream & out, std::string & buff_i)
{
        std::stringstream ss;
        ss << out.rdbuf();
        std::string myString = ss.str();
        return buff_i;
}
