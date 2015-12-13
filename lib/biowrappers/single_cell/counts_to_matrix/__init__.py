import tasks

def convert_counts_to_matrix_pipeline(sch, in_file, out_file, col_id_fields, data_type, filters, row_id_fields):
    sch.transform('convert_counts_to_matrix', in_file.axes, {},
                  tasks.convert_counts_to_input_matrix,
                  None,
                  in_file.input,
                  out_file.output,
                  col_id_fields,
                  data_type,
                  filters,
                  row_id_fields)