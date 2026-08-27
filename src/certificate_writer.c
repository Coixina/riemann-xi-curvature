#include "certificate_writer.h"

#include <limits.h>


void
certificate_writer_init(
    certificate_writer_t *writer
)
{
    writer->file = NULL;
    writer->leaf_count_offset = -1;
    writer->leaf_count = 0;
    writer->is_open = 0;
    writer->is_finalized = 0;
}


int
certificate_writer_open(
    certificate_writer_t *writer,
    const char *path,
    ulong large_t_N,
    const fmpz_t kappa_num,
    const fmpz_t kappa_den,
    slong generator_precision
)
{
    FILE *file;

    long count_offset;


    if (writer == NULL ||
        path == NULL)
    {
        return 0;
    }

    if (writer->is_open)
        return 0;

    if (large_t_N < 2)
        return 0;

    if (fmpz_sgn(kappa_den) <= 0)
        return 0;

    if (generator_precision <= 0)
        return 0;


    file =
        fopen(
            path,
            "w+"
        );

    if (file == NULL)
        return 0;


    if (fprintf(
            file,
            "RH_FINITE_CERTIFICATE %d\n",
            RH_CERTIFICATE_FORMAT_VERSION
        ) < 0)
    {
        fclose(file);
        return 0;
    }


    if (fprintf(
            file,
            "large_t_N %lu\n",
            large_t_N
        ) < 0)
    {
        fclose(file);
        return 0;
    }


    if (fprintf(
            file,
            "kappa_num "
        ) < 0)
    {
        fclose(file);
        return 0;
    }

    if (fmpz_fprint(
            file,
            kappa_num
        ) <= 0)
    {
        fclose(file);
        return 0;
    }

    if (fputc('\n', file) == EOF)
    {
        fclose(file);
        return 0;
    }


    if (fprintf(
            file,
            "kappa_den "
        ) < 0)
    {
        fclose(file);
        return 0;
    }

    if (fmpz_fprint(
            file,
            kappa_den
        ) <= 0)
    {
        fclose(file);
        return 0;
    }

    if (fputc('\n', file) == EOF)
    {
        fclose(file);
        return 0;
    }


    if (fprintf(
            file,
            "generator_precision %ld\n",
            generator_precision
        ) < 0)
    {
        fclose(file);
        return 0;
    }


    if (fprintf(
            file,
            "leaf_count "
        ) < 0)
    {
        fclose(file);
        return 0;
    }


    count_offset =
        ftell(file);

    if (count_offset < 0)
    {
        fclose(file);
        return 0;
    }


    if (fprintf(
            file,
            "%020lu\n",
            0UL
        ) < 0)
    {
        fclose(file);
        return 0;
    }


    if (fprintf(
            file,
            "BEGIN_LEAVES\n"
        ) < 0)
    {
        fclose(file);
        return 0;
    }


    writer->file =
        file;

    writer->leaf_count_offset =
        count_offset;

    writer->leaf_count = 0;

    writer->is_open = 1;
    writer->is_finalized = 0;


    return 1;
}


int
certificate_writer_write_leaf(
    certificate_writer_t *writer,
    const interval_leaf_t *leaf
)
{
    if (writer == NULL ||
        leaf == NULL)
    {
        return 0;
    }

    if (!writer->is_open ||
        writer->is_finalized ||
        writer->file == NULL)
    {
        return 0;
    }

    if (!interval_leaf_is_valid(leaf))
        return 0;

    if (writer->leaf_count == ULONG_MAX)
        return 0;


    if (fprintf(
            writer->file,
            "%lu ",
            leaf->depth
        ) < 0)
    {
        return 0;
    }


    if (fmpz_fprint(
            writer->file,
            leaf->index
        ) <= 0)
    {
        return 0;
    }


    if (fputc(
            '\n',
            writer->file
        ) == EOF)
    {
        return 0;
    }


    writer->leaf_count++;

    return 1;
}


int
certificate_writer_finalize(
    certificate_writer_t *writer
)
{
    long end_position;


    if (writer == NULL)
        return 0;

    if (!writer->is_open ||
        writer->is_finalized ||
        writer->file == NULL)
    {
        return 0;
    }


    if (fprintf(
            writer->file,
            "END_LEAVES\n"
        ) < 0)
    {
        return 0;
    }


    if (fflush(
            writer->file
        ) != 0)
    {
        return 0;
    }


    end_position =
        ftell(
            writer->file
        );

    if (end_position < 0)
        return 0;


    if (fseek(
            writer->file,
            writer->leaf_count_offset,
            SEEK_SET
        ) != 0)
    {
        return 0;
    }


    if (fprintf(
            writer->file,
            "%020lu",
            writer->leaf_count
        ) < 0)
    {
        return 0;
    }


    if (fflush(
            writer->file
        ) != 0)
    {
        return 0;
    }


    if (fseek(
            writer->file,
            end_position,
            SEEK_SET
        ) != 0)
    {
        return 0;
    }


    if (fclose(
            writer->file
        ) != 0)
    {
        writer->file = NULL;
        writer->is_open = 0;
        return 0;
    }


    writer->file = NULL;
    writer->is_open = 0;
    writer->is_finalized = 1;

    return 1;
}


void
certificate_writer_abort(
    certificate_writer_t *writer
)
{
    if (writer == NULL)
        return;

    if (writer->file != NULL)
    {
        fclose(
            writer->file
        );
    }

    writer->file = NULL;
    writer->leaf_count_offset = -1;
    writer->is_open = 0;
    writer->is_finalized = 0;
}