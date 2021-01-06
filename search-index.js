var searchIndex = JSON.parse('{\
"kff":{"doc":"kff-rust-apiA rust API to read and write Kmer File Format…","i":[[0,"data","kff","Declaration of trait to read and write kmer section",null,null],[8,"Reader","kff::data","Reader trait must be implement by Struct want read kmer…",null,null],[10,"k","","Get size of kmers",0,[[]]],[10,"input","","Get a reference to the read stream",0,[[]]],[10,"max_kmer","","Get the max number of kmer in a block",0,[[]]],[10,"block_n","","Get the actual number of kmer remains in block",0,[[]]],[10,"data_size","","Get the size of data attach to each kmer in bytes",0,[[]]],[10,"read_block","","Read on kmer block and return number of bytes read",0,[[],["result",6]]],[10,"block_seq","","Get the sequence of actual block",0,[[],["seq2slice",6]]],[10,"block_data","","Get the data of actutal block",0,[[]]],[10,"rev_encoding","","Get bit used to perform reverse encoding",0,[[]]],[10,"decrease_n","","Reduce block_n by one",0,[[]]],[11,"read_n","","Read the number of kmer in block",0,[[],["result",6]]],[11,"read_seq","","Read sequence of actual block need number of nucleotide…",0,[[],[["result",6],["bitbox",6]]]],[11,"read_data","","Read data of actual block",0,[[],[["vec",3],["result",6]]]],[11,"next_kmer","","Get the next kmer",0,[[],[["result",6],["kmer",3]]]],[8,"Writer","","Write trait must be implement by Struct want write kmer…",null,null],[10,"data_size","","Get the size of data attach to each kmer in bytes",1,[[]]],[10,"max","","Get the max number of kmer in a block",1,[[]]],[10,"nb_block_offset","","Position of nb_block_offset need to write the number of…",1,[[]]],[10,"nb_block","","Number of block in this section",1,[[]]],[10,"is_close","","Return True if section is close don\'t write any other…",1,[[]]],[10,"output","","Get a reference of output",1,[[]]],[10,"set_nb_block","","Change number of block in section",1,[[]]],[10,"set_close","","Change state of [Writer::is_close]",1,[[]]],[10,"nb_kmer","","Compute the number of kmer from number of sequence",1,[[]]],[11,"close","","Close the section, can return an error because number of…",1,[[],["result",6]]],[11,"increment_nb_block","","Increase the number of block",1,[[],["result",6]]],[11,"check_block","","Verify the number of the number of kmer and data length…",1,[[],["result",6]]],[0,"error","kff","Declaration of error type",null,null],[4,"Header","kff::error","Error can be generate in Header section",null,null],[13,"BadEncoding","","Value associate to each nucleotide must be differente",2,null],[4,"Variables","","Error can be generate in Variable section",null,null],[13,"KMissing","","To read Raw and Minimizer section, variable k must be set",3,null],[13,"MMissing","","To read Minimizer section, variable m must be set",3,null],[13,"MaxMissing","","To read Raw and Minimizer section, variable max must be set",3,null],[13,"DataSizeMissing","","To read Raw and Minimizer section variable data_size must…",3,null],[4,"Data","","Error can be generate in Raw or Minimizer section",null,null],[13,"NUpperThanMax","","Number of kmer in block can\'t be larger than max",4,null],[13,"NbKmerNbDataDiff","","Number of kmer and number of data must be equal",4,null],[13,"ToManyBlock","","We can\'t have more than u32::max() block in section",4,null],[4,"Minimizer","","Error can be generate in Minimizer section",null,null],[13,"MinimizerSizeMDiff","","Size of minimizer sequence is differente than m variables",5,null],[4,"Kff","","Error can be generate at any moment durring kff reading or…",null,null],[13,"UnknowSectionType","","This api support only Variables, Raw and Minimizer section",6,null],[13,"NotSupportVersionNumber","","Not suport version of KFF file",6,null],[4,"Error","","Generale error type use in this parser",null,null],[13,"Header","","",7,null],[13,"Variables","","",7,null],[13,"Data","","",7,null],[13,"Minimizer","","",7,null],[13,"Kff","","",7,null],[0,"kff","kff","Declaration of Kff::Reader and Kff::Writer",null,null],[3,"Reader","kff::kff","A Kff Reader",null,null],[3,"Writer","","A Kff Writer",null,null],[11,"new","","Create a new Kff reader from a reading stream",8,[[],["result",6]]],[11,"major","","Get the major verion number",8,[[]]],[11,"minor","","Get the minor verion number",8,[[]]],[11,"encoding","","Get encoding used",8,[[]]],[11,"rev_encoding","","Get reverse encoding used",8,[[]]],[11,"comment","","Get comment",8,[[]]],[11,"input","","Get a mutable reference to input stream",8,[[]]],[11,"variables","","Get a mutable reference of global variable",8,[[],["variables",6]]],[11,"next_section","","Get the next kmer section we have to parse",8,[[],[["box",3],["result",6]]]],[11,"new","","Create a Kff Writer from a seekable output, encoding and…",9,[[],["result",6]]],[11,"encoding","","Get encoding used",9,[[]]],[11,"variables","","Get a mutable reference of global variable",9,[[],["variables",6]]],[11,"write_variables","","Write variable section",9,[[],["result",6]]],[11,"write_raw_section","","Write a raw section, with sequence encode in 2 bits",9,[[],["result",6]]],[11,"write_raw_seq_section","","Write a raw section with sequence encode in ASCII",9,[[],["result",6]]],[11,"write_minimizer_section","","Write a minimizer section, with sequence encode in 2 bits",9,[[],["result",6]]],[11,"write_minimizer_seq_section","","Write a raw section with sequence not encode in ASCII",9,[[],["result",6]]],[0,"kmer","kff","Define a kmer type",null,null],[3,"Kmer","kff::kmer","Kmer type",null,null],[11,"new","","Create a new kmer",10,[[["box",3],["seq2bits",6]]]],[11,"seq","","Get sequence in 2 bits encoding",10,[[],["seq2bits",6]]],[11,"data","","Get data",10,[[]]],[11,"len","","Get length of kmer",10,[[]]],[11,"is_empty","","Return true if kmer length is 0",10,[[]]],[0,"minimizer","kff","Declaration of Minimizer section Reader and Writer",null,null],[3,"Reader","kff::minimizer","A Minimizer section reader implement [data::Reader]",null,null],[3,"Writer","","A Minimizer section writer implement [data::Writer]",null,null],[11,"new","","Create a new reader with a reference of kff::Reader",11,[[["kffreader",3]],["result",6]]],[11,"new","","Create a new Minimizer section writer",12,[[["variables",6]],["result",6]]],[11,"write_block","","Write a minimizer block",12,[[["bitslice",6]],["result",6]]],[11,"write_seq_block","","Write a minimizer block, where sequence is encode in ASCII",12,[[],["result",6]]],[0,"raw","kff","Declaration of Raw section Reader and Writer",null,null],[3,"Reader","kff::raw","A Raw section reader implement [data::Reader]",null,null],[3,"Writer","","A Raw section writer implement [data::Writer]",null,null],[11,"new","","Create a new reader with a reference of kff::Reader",13,[[["kffreader",3]],["result",6]]],[11,"new","","Create a new Raw section writer",14,[[["variables",6]],["result",6]]],[11,"write_block","","Write a raw block",14,[[["bitslice",6]],["result",6]]],[11,"write_seq_block","","Write a raw block, where sequence is encode in ASCII",14,[[],["result",6]]],[0,"seq2bits","kff","Declaration of type and trait to store 2 bits…",null,null],[6,"Seq2Slice","kff::seq2bits","Syntaxic sugar",null,null],[6,"Seq2Bits","","Syntaxic sugar",null,null],[8,"Nuc2Bits","","Build Self, from some source",null,null],[10,"from_nuc","","Build Self from nucleotide encode in ASCII",15,[[]]],[10,"from_bitslice","","Build Self from nucleotide encode on two bits store in…",15,[[["bitslice",6]]]],[10,"from_bits","","Build Self from nucleotide encode on two bits store in…",15,[[]]],[8,"Bits2Nuc","","Convert Self into ASCII encoding",null,null],[10,"into_nuc","","Convert Self in ASCII encoding",16,[[],["box",3]]],[8,"RangeNuc","","",null,null],[10,"range_nuc","","",17,[[["range",3]],["seq2slice",6]]],[0,"utils","kff","",null,null],[6,"Order","kff::utils","Read order of bytes in file",null,null],[6,"BitOrd","","Order of bit for bitvec",null,null],[6,"BitBox","","Syntaxic sugar around bitvec::BitBox",null,null],[6,"BitVec","","Syntaxic sugar around bitvec::BitVec",null,null],[6,"BitSlice","","Syntaxic sugar around bitvec::BitSlice",null,null],[0,"variables","kff","Management of global variable",null,null],[6,"Variables","kff::variables","Variables is a specialisation of HashMap",null,null],[8,"Reader","","Trait to read global variable section",null,null],[10,"deserialize","","Read variable from input",18,[[],["result",6]]],[8,"Writer","","Trait to write global variable section",null,null],[10,"serialize","","Write variable in output",19,[[],["result",6]]],[8,"Variables1","","Trait of variable needed by KFF 1.0",null,null],[10,"k","","Get value of k",20,[[],["result",6]]],[10,"m","","Get value of m",20,[[],["result",6]]],[10,"max","","Get value of max",20,[[],["result",6]]],[10,"data_size","","Get value of data_size",20,[[],["result",6]]],[11,"from","kff::error","",2,[[]]],[11,"into","","",2,[[]]],[11,"to_string","","",2,[[],["string",3]]],[11,"borrow","","",2,[[]]],[11,"borrow_mut","","",2,[[]]],[11,"try_from","","",2,[[],["result",4]]],[11,"try_into","","",2,[[],["result",4]]],[11,"type_id","","",2,[[],["typeid",3]]],[11,"from","","",3,[[]]],[11,"into","","",3,[[]]],[11,"to_string","","",3,[[],["string",3]]],[11,"borrow","","",3,[[]]],[11,"borrow_mut","","",3,[[]]],[11,"try_from","","",3,[[],["result",4]]],[11,"try_into","","",3,[[],["result",4]]],[11,"type_id","","",3,[[],["typeid",3]]],[11,"from","","",4,[[]]],[11,"into","","",4,[[]]],[11,"to_string","","",4,[[],["string",3]]],[11,"borrow","","",4,[[]]],[11,"borrow_mut","","",4,[[]]],[11,"try_from","","",4,[[],["result",4]]],[11,"try_into","","",4,[[],["result",4]]],[11,"type_id","","",4,[[],["typeid",3]]],[11,"from","","",5,[[]]],[11,"into","","",5,[[]]],[11,"to_string","","",5,[[],["string",3]]],[11,"borrow","","",5,[[]]],[11,"borrow_mut","","",5,[[]]],[11,"try_from","","",5,[[],["result",4]]],[11,"try_into","","",5,[[],["result",4]]],[11,"type_id","","",5,[[],["typeid",3]]],[11,"from","","",6,[[]]],[11,"into","","",6,[[]]],[11,"to_string","","",6,[[],["string",3]]],[11,"borrow","","",6,[[]]],[11,"borrow_mut","","",6,[[]]],[11,"try_from","","",6,[[],["result",4]]],[11,"try_into","","",6,[[],["result",4]]],[11,"type_id","","",6,[[],["typeid",3]]],[11,"from","","",7,[[]]],[11,"into","","",7,[[]]],[11,"to_string","","",7,[[],["string",3]]],[11,"borrow","","",7,[[]]],[11,"borrow_mut","","",7,[[]]],[11,"try_from","","",7,[[],["result",4]]],[11,"try_into","","",7,[[],["result",4]]],[11,"type_id","","",7,[[],["typeid",3]]],[11,"from","kff::kff","",8,[[]]],[11,"into","","",8,[[]]],[11,"borrow","","",8,[[]]],[11,"borrow_mut","","",8,[[]]],[11,"try_from","","",8,[[],["result",4]]],[11,"try_into","","",8,[[],["result",4]]],[11,"type_id","","",8,[[],["typeid",3]]],[11,"from","","",9,[[]]],[11,"into","","",9,[[]]],[11,"borrow","","",9,[[]]],[11,"borrow_mut","","",9,[[]]],[11,"try_from","","",9,[[],["result",4]]],[11,"try_into","","",9,[[],["result",4]]],[11,"type_id","","",9,[[],["typeid",3]]],[11,"from","kff::kmer","",10,[[]]],[11,"into","","",10,[[]]],[11,"borrow","","",10,[[]]],[11,"borrow_mut","","",10,[[]]],[11,"try_from","","",10,[[],["result",4]]],[11,"try_into","","",10,[[],["result",4]]],[11,"type_id","","",10,[[],["typeid",3]]],[11,"from","kff::minimizer","",11,[[]]],[11,"into","","",11,[[]]],[11,"into_iter","","",11,[[]]],[11,"borrow","","",11,[[]]],[11,"borrow_mut","","",11,[[]]],[11,"try_from","","",11,[[],["result",4]]],[11,"try_into","","",11,[[],["result",4]]],[11,"type_id","","",11,[[],["typeid",3]]],[11,"from","","",12,[[]]],[11,"into","","",12,[[]]],[11,"borrow","","",12,[[]]],[11,"borrow_mut","","",12,[[]]],[11,"try_from","","",12,[[],["result",4]]],[11,"try_into","","",12,[[],["result",4]]],[11,"type_id","","",12,[[],["typeid",3]]],[11,"from","kff::raw","",13,[[]]],[11,"into","","",13,[[]]],[11,"into_iter","","",13,[[]]],[11,"borrow","","",13,[[]]],[11,"borrow_mut","","",13,[[]]],[11,"try_from","","",13,[[],["result",4]]],[11,"try_into","","",13,[[],["result",4]]],[11,"type_id","","",13,[[],["typeid",3]]],[11,"from","","",14,[[]]],[11,"into","","",14,[[]]],[11,"borrow","","",14,[[]]],[11,"borrow_mut","","",14,[[]]],[11,"try_from","","",14,[[],["result",4]]],[11,"try_into","","",14,[[],["result",4]]],[11,"type_id","","",14,[[],["typeid",3]]],[11,"k","kff::minimizer","",11,[[]]],[11,"input","","",11,[[]]],[11,"max_kmer","","",11,[[]]],[11,"block_n","","",11,[[]]],[11,"data_size","","",11,[[]]],[11,"block_seq","","",11,[[],["bitslice",6]]],[11,"block_data","","",11,[[]]],[11,"rev_encoding","","",11,[[]]],[11,"decrease_n","","",11,[[]]],[11,"read_block","","",11,[[],["result",6]]],[11,"k","kff::raw","",13,[[]]],[11,"input","","",13,[[]]],[11,"max_kmer","","",13,[[]]],[11,"block_n","","",13,[[]]],[11,"data_size","","",13,[[]]],[11,"block_seq","","",13,[[],["bitslice",6]]],[11,"block_data","","",13,[[]]],[11,"rev_encoding","","",13,[[]]],[11,"decrease_n","","",13,[[]]],[11,"read_block","","",13,[[],["result",6]]],[11,"data_size","kff::minimizer","",12,[[]]],[11,"max","","",12,[[]]],[11,"nb_block_offset","","",12,[[]]],[11,"nb_block","","",12,[[]]],[11,"is_close","","",12,[[]]],[11,"output","","",12,[[]]],[11,"set_nb_block","","",12,[[]]],[11,"set_close","","",12,[[]]],[11,"nb_kmer","","",12,[[]]],[11,"data_size","kff::raw","",14,[[]]],[11,"max","","",14,[[]]],[11,"nb_block_offset","","",14,[[]]],[11,"nb_block","","",14,[[]]],[11,"is_close","","",14,[[]]],[11,"output","","",14,[[]]],[11,"set_nb_block","","",14,[[]]],[11,"set_close","","",14,[[]]],[11,"nb_kmer","","",14,[[]]],[11,"from_nuc","kff","",21,[[]]],[11,"from_bitslice","","",21,[[["bitslice",6]]]],[11,"from_bits","","",21,[[]]],[11,"into_nuc","","",21,[[],["box",3]]],[11,"into_nuc","","",22,[[],["box",3]]],[11,"range_nuc","","",21,[[["range",3]],["seq2slice",6]]],[11,"deserialize","","",23,[[],["result",6]]],[11,"serialize","","",23,[[],["result",6]]],[11,"k","","",23,[[],["result",6]]],[11,"m","","",23,[[],["result",6]]],[11,"max","","",23,[[],["result",6]]],[11,"data_size","","",23,[[],["result",6]]],[11,"drop","kff::minimizer","",12,[[]]],[11,"drop","kff::raw","",14,[[]]],[11,"next","kff::minimizer","",11,[[],["option",4]]],[11,"next","kff::raw","",13,[[],["option",4]]],[11,"eq","kff::kmer","",10,[[["kmer",3]]]],[11,"ne","","",10,[[["kmer",3]]]],[11,"fmt","kff::error","",2,[[["formatter",3]],["result",6]]],[11,"fmt","","",3,[[["formatter",3]],["result",6]]],[11,"fmt","","",4,[[["formatter",3]],["result",6]]],[11,"fmt","","",5,[[["formatter",3]],["result",6]]],[11,"fmt","","",6,[[["formatter",3]],["result",6]]],[11,"fmt","","",7,[[["formatter",3]],["result",6]]],[11,"fmt","kff::kmer","",10,[[["formatter",3]],["result",6]]],[11,"fmt","kff::error","",2,[[["formatter",3]],["result",6]]],[11,"fmt","","",3,[[["formatter",3]],["result",6]]],[11,"fmt","","",4,[[["formatter",3]],["result",6]]],[11,"fmt","","",5,[[["formatter",3]],["result",6]]],[11,"fmt","","",6,[[["formatter",3]],["result",6]]],[11,"fmt","","",7,[[["formatter",3]],["result",6]]],[11,"source","","",7,[[],[["option",4],["error",8]]]]],"p":[[8,"Reader"],[8,"Writer"],[4,"Header"],[4,"Variables"],[4,"Data"],[4,"Minimizer"],[4,"Kff"],[4,"Error"],[3,"Reader"],[3,"Writer"],[3,"Kmer"],[3,"Reader"],[3,"Writer"],[3,"Reader"],[3,"Writer"],[8,"Nuc2Bits"],[8,"Bits2Nuc"],[8,"RangeNuc"],[8,"Reader"],[8,"Writer"],[8,"Variables1"],[6,"Seq2Bits"],[6,"Seq2Slice"],[6,"Variables"]]}\
}');
addSearchOptions(searchIndex);initSearch(searchIndex);