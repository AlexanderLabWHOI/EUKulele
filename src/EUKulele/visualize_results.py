'''
Create plots to convey taxonomy.
'''

import os
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

#.loc[[name == curr for name in final_frame.loc[name_level]],["Sum"]]

def countClassifs(level, level_hierarchy, name_level, df):
    ''' Count the classifications assigned to each level and type. '''

    set_list = set()

    classifications = list(df.loc[df["classification_level"] == level]["classification"])
    counts = list(df.loc[df["classification_level"] == level]["counts"])
    transcript_names = list(df.loc[df["classification_level"] == level]["transcript_name"])
    print(level_hierarchy, flush=True)
    match_loc = int(np.where([curr == level for curr in level_hierarchy])[0])
    print(match_loc, flush=True)
    for curr in range(match_loc + 1,len(level_hierarchy)):
        classification_curr = list(df.loc[df["classification_level"] == \
                                          level_hierarchy[curr]]["full_classification"])
        transcripts_curr = list(df.loc[df["classification_level"] == \
                                       level_hierarchy[curr]]["transcript_name"])
        correct_index = list(np.where([len(str(cr).split(";")) >= \
                                       abs(-1-(curr-match_loc)) \
                                       for cr in classification_curr])[0])

        classification_curr = [classification_curr[cr2] for cr2 in correct_index]
        classifs_curr = [str(cr).split(";")[match_loc].strip() \
                         for cr in classification_curr]
        transcripts_curr = [transcripts_curr[cr2] for cr2 in correct_index]
        counts_curr = list(df.loc[df["classification_level"] == \
                                  level_hierarchy[curr]]["counts"])
        counts_curr = [counts_curr[cr2] for cr2 in correct_index]

        # add to running list
        classifications.extend(classifs_curr)
        transcript_names.extend(transcripts_curr)
        counts.extend(counts_curr)
    lostducklings = df[-df['transcript_name'].isin(transcript_names)]
    counts.extend(list(lostducklings.counts))
    transcript_names.extend(list(lostducklings.transcript_name))
    classifications.extend(["NoClassification"] * len(lostducklings.index))

    classifications = [str(cr).strip().strip(""''"]['") for cr in classifications]
    full_list = classifications
    set_list.update(set(classifications))

    transcripts_classes = pd.DataFrame({"classifications": classifications,\
                                        "transcript_names": transcript_names})
    transcripts_classes = transcripts_classes.groupby("classifications").agg(\
        {"transcript_names": lambda x: ';'.join(x)})

    # Apparently now when you groupby, the grouping becomes the index.
    if transcripts_classes.index.name != "classifications":
        transcripts_classes = transcripts_classes.set_index('classifications')
    transcripts_classes = transcripts_classes.loc[sorted(list(set_list))]

    full_frame = pd.DataFrame({name_level: classifications, "Counts": counts})
    final_frame = full_frame.groupby(name_level)['Counts'].agg(Sum='sum', Count='count')
    # gives both the sum of cts & # transcripts
    final_frame["GroupedTranscripts"] = list(transcripts_classes.transcript_names)

    end_frame = pd.DataFrame({name_level: sorted(list(set_list)),
                                "Counts": tuple([float(final_frame[final_frame.index == \
                                                                   curr].Sum) \
                                           for curr in sorted(list(set_list))]),
                                "NumTranscripts": [full_list.count(curr) for curr \
                                                   in sorted(list(set_list))],
                                "GroupedTranscripts": list(transcripts_classes.transcript_names)})
    return classifications, end_frame

def countClassifsNoCounts(level, level_hierarchy, name_level, df):
    ''' Same counting procedure, but sans Salmon counts. '''

    set_list = set()

    classifications = list(df.loc[df["classification_level"] == level]["classification"])
    transcript_names = list(df.loc[df["classification_level"] == level]["transcript_name"])
    print(level_hierarchy, flush=True)
    print(level, flush=True)
    if len(np.where([curr == level.lower() for curr in level_hierarchy])) == 0:
        return None, None
    match_loc = int(np.where([curr == level.lower() for curr in level_hierarchy])[0])
    print(match_loc, flush=True)

    for curr in range(match_loc + 1,len(level_hierarchy)):
        ## MAKE THE TWO CHANGES FROM ABOVE HERE!!
        classification_curr = list(df.loc[df["classification_level"] == \
                                          level_hierarchy[curr]]["full_classification"])
        transcripts_curr = list(df.loc[df["classification_level"] == \
                                       level_hierarchy[curr]]["transcript_name"])
        correct_index = list(np.where([len(str(cr).split(";")) >= \
                                       abs(match_loc) \
                                       #abs(-1-(curr-match_loc)) \
                                       for cr in classification_curr])[0])

        classification_curr = [classification_curr[cr2] for cr2 in correct_index]
        transcripts_curr = [transcripts_curr[cr2] for cr2 in correct_index]
        classifs_curr = [str(cr).split(";")[match_loc].strip() \
                         for cr in classification_curr]
        #classifs_curr = [str(cr).split(";")[-1-(curr-match_loc)].strip() \
        #                 for cr in classification_curr]

        transcript_names.extend(transcripts_curr)
        classifications.extend(classifs_curr)

    # a vector containing all of the classifications given at the current taxonomic level
    classifications = [str(cr).strip().strip(""''"]['") for cr in classifications]
    full_list = classifications
    set_list.update(set(classifications))

    transcripts_classes = pd.DataFrame({"classifications": classifications,
                                        "transcript_names": transcript_names})
    transcripts_classes = transcripts_classes.groupby("classifications").\
                            agg({"transcript_names": lambda x: ';'.join(x)})
    transcripts_classes.sort_values(by = "classifications", inplace=True)

    final_frame = pd.DataFrame({name_level: sorted(list(set_list)),
                                "Counts": [full_list.count(curr) for \
                                           curr in sorted(list(set_list))],
                                "GroupedTranscripts": list(transcripts_classes.\
                                                           transcript_names)})

    return classifications, final_frame

def stripClassifData(df, use_counts, level_hierarchy):
    ''' Pull the classification data from the tokenized taxonomy file. '''

    #level_hierarchy = ['supergroup','division','class','order',\
    #                   'family','genus','species']

    return_dict_list = dict()
    return_dict_frame = dict()
    for curr_level in level_hierarchy:
        if use_counts == True:
            curr_list, curr_counts = countClassifs(curr_level, level_hierarchy,
                                                   curr_level.capitalize(), df)
            return_dict_list[curr_level] = curr_list
            return_dict_frame[curr_level] = curr_counts
        else:
            curr_list, curr_counts = countClassifsNoCounts(curr_level,\
                                                           level_hierarchy,\
                                                           curr_level.capitalize(),\
                                                           df)
            return_dict_list[curr_level] = curr_list
            return_dict_frame[curr_level] = curr_counts
    return return_dict_list, return_dict_frame

def createPlotDataFrame(curr_df_start, cutoff_relative = 0.1,
                        transcript_or_counts="NumTranscripts"):
    ''' Creates a consolidated dataframe for use in plotting. '''

    ## CREATE AGGREGATED COUNTS BY SAMPLE ##
    curr_df_summed = curr_df_start.groupby("Sample")[transcript_or_counts].agg(AllCts='sum')
    curr_df_summed = curr_df_start.join(curr_df_summed,how="left",on="Sample")

    ## CALCULATE RELATIVE COUNTS ##
    curr_df_summed["Rel_Counts"] = curr_df_summed[transcript_or_counts] / curr_df_summed["AllCts"]
    curr_df_plot = curr_df_summed[["OfInterest","Sample","Rel_Counts"]]
    pivoted = curr_df_plot.pivot(index='Sample', columns='OfInterest', values='Rel_Counts')
    pivoted = pivoted.reindex(sorted(pivoted.columns), axis=1)

    ## If we're plotting both counts and transcripts, still decide on cutoff via transcripts ##
    if transcript_or_counts == "Counts":
        curr_df_summed = curr_df_start.groupby("Sample")["NumTranscripts"].agg(AllCts='sum')
        curr_df_summed = curr_df_start.join(curr_df_summed,how="left",on="Sample")
        curr_df_summed["Rel_Counts_Transcripts"] = curr_df_summed["NumTranscripts"]/ \
                                                   curr_df_summed["AllCts"]
        curr_df_summed = curr_df_summed[["OfInterest","Sample","Rel_Counts_Transcripts"]]
        pivoted_on_transcripts = curr_df_summed.pivot(index='Sample',
                                                      columns='OfInterest',
                                                      values='Rel_Counts_Transcripts')
        pivoted_on_transcripts = pivoted_on_transcripts.reindex(sorted(pivoted_on_transcripts.\
                                                                       columns), axis=1)
        pivoted_agg = list(pivoted_on_transcripts.max(axis = 0, skipna = True))
    else:
        pivoted_agg = list(pivoted.max(axis = 0, skipna = True))

    ## Leave the columns that meet the threshold; sum the others into an "Other" column ##
    chosen_cols = [curr for curr in range(len(pivoted_agg)) if pivoted_agg[curr] > \
                   cutoff_relative]
    chosen_cols_other = [curr for curr in range(len(pivoted_agg)) if pivoted_agg[curr] \
                         <= cutoff_relative]
    pivoted_agg = list(pivoted.iloc[:,chosen_cols_other].sum(axis = 1,\
                                                             skipna = True))

    ## Modify the output dataframe accordingly ##
    pivoted = pivoted.iloc[:,chosen_cols]
    pivoted["Other"] = pivoted_agg

    return pivoted

def makeConcatFrame(curr_df, new_df, level, sample_name, use_counts):
    ''' Add current DataFrame to running list. '''

    new_df = pd.DataFrame(new_df)
    if new_df.empty:
        return curr_df
    if use_counts == True:
        new_df = pd.DataFrame(new_df.reset_index())
        new_df = new_df.drop(columns = "index")
        new_df.columns = [level,"Counts","NumTranscripts","GroupedTranscripts"]
    else:
        new_df = pd.DataFrame(new_df)
        new_df.columns = [level,"NumTranscripts","GroupedTranscripts"]

    new_df["Sample"] = sample_name
    return pd.concat([curr_df, new_df], sort=True)

def visualize_all_results(out_prefix, out_dir, est_dir, samples_dir,
                          prot_extension, nucle_extension, use_counts,
                          rerun, level_hierarchy, core = False):
    ''' Main visualization method. '''

    results_frame = dict()
    results_counts_dir = os.path.join(out_dir, "taxonomy_counts")
    results_viz_dir = os.path.join(out_dir, "taxonomy_visualization")
    if core:
        results_counts_dir = os.path.join(out_dir, "core_taxonomy_counts")
        results_viz_dir = os.path.join(out_dir, "core_taxonomy_visualization")

    ### READ IN RESULTS FILES FROM MET DIR THAT FIT SAMPLE SPEC FROM CONFIG ###
    samples = os.listdir(samples_dir)
    good_samples = 0
    for s_curr in samples:
        file_name = ".".join(s_curr.split(".")[0:-1]) + "-estimated-taxonomy.out"
        if (prot_extension in s_curr.split(".")[-1]) | \
           (nucle_extension in s_curr.split(".")[-1]):
            if not os.path.isfile(os.path.join(est_dir, file_name)):
                print("One of the files, " + s_curr + \
                      ", in the sample directory did not complete successfully.")
                if not core:
                    sys.exit(1)
            else:
                results_frame[file_name] = pd.read_csv(os.path.join(est_dir,\
                                                                    file_name),\
                                                       sep = "\t", index_col=0)
        good_samples = good_samples + 1

    if (good_samples == 0) & (not core):
        print("No taxonomic estimation files found! Exiting...")
        sys.exit(1)

    list_results = dict()
    frame_results = dict()

    # characterizing by major classes
    for curr in results_frame.keys():
        list_results[curr], frame_results[curr] = stripClassifData(results_frame[curr],
                                                                   use_counts,
                                                                   level_hierarchy)

    counts_all = dict()
    #level_hierarchy = ['supergroup','division','class','order','family',
    #                   'genus','species']
    for l_curr in level_hierarchy:
        counts_all[l_curr] = pd.DataFrame(columns = [l_curr.capitalize(),
                                                     "NumTranscripts",
                                                     "GroupedTranscripts",
                                                     "Sample"])

    for curr in results_frame.keys():
        if results_frame[curr] is None:
            next
        sample_name = curr #curr.split("-")[0]
        for l in level_hierarchy:
            curr_df = counts_all[l]
            new_df = frame_results[curr][l]
            counts_all[l] = makeConcatFrame(curr_df, new_df, l.capitalize(),
                                            sample_name, use_counts)

    for l in level_hierarchy:
        ### SAVE THE CSVs OF THE DATA ###
        prefix = out_prefix
        os.system("mkdir -p " + results_counts_dir)
        counts_all[l].to_csv(os.path.join(results_counts_dir, prefix +\
                                          "_all_" + l + "_counts.csv"))

        if not os.path.isfile(os.path.join(results_counts_dir, prefix +\
                                            "_all_" + l + "_counts.csv")):
            print("Taxonomy counts were not successfully generated. Check log for details.")
            sys.exit(1)

        ### INITIALIZE VARIABLES FOR LOOP ###
        Curr_Variable = l.capitalize()
        cutoff_counts = 100000
        cutoff_transcripts = 100

        ### GRAB DATAFRAME FROM SAVED DATA ###
        curr_df_start = counts_all[l].reset_index()
        if Curr_Variable in curr_df_start:
            curr_df_start["OfInterest"] = curr_df_start[Curr_Variable]
        else:
            curr_df_start["OfInterest"] = []

        sns.set()

        ## CREATE AGGREGATED COUNTS BY SAMPLE ##
        curr_df_summed = curr_df_start.groupby("Sample")["NumTranscripts"].agg(AllCts='sum')
        curr_df_summed = curr_df_start.join(curr_df_summed,how="left",on="Sample")

        ## CALCULATE RELATIVE COUNTS ##
        curr_df_summed["Rel_Counts"] = curr_df_summed["NumTranscripts"] / curr_df_summed["AllCts"]
        curr_df_plot = curr_df_summed[["OfInterest","Sample","Rel_Counts"]]
        pivoted = curr_df_plot.pivot(index='Sample', columns='OfInterest', values='Rel_Counts')

        pivoted_agg = list(pivoted.sum(axis = 0, skipna = True))
        chosen_cols = [curr for curr in range(len(pivoted_agg)) if pivoted_agg[curr] > 0.1]
        pivoted = pivoted.iloc[:,chosen_cols]

        c25 = ["dodgerblue2", "#E31A1C",\
          "green4",\
          "#6A3D9A",\
          "#FF7F00",\
          "black", "gold1",\
          "skyblue2", "#FB9A99",\
          "palegreen2",\
          "#CAB2D6",\
          "#FDBF6F",\
          "gray70", "khaki2",\
          "maroon", "orchid1", "deeppink1", "blue1", "steelblue4",\
          "darkturquoise", "green1", "yellow4", "yellow3",\
          "darkorange4", "brown"]

        fig = plt.figure(figsize=(15,7.5))
        if len(set(curr_df_start["OfInterest"])) == 0:
            continue
        sns_palette = sns.palplot(sns.color_palette("Set1",\
                                                    n_colors=len(set(curr_df_start["OfInterest"]))))

        ### CREATE PLOTS ###
        if use_counts == False:
            fig = plt.figure(figsize=(15,7.5))
            fig.set_facecolor('white')
            ax = fig.add_axes([0.1, 0.1, 0.6, 0.75])
            pivoted = createPlotDataFrame(curr_df_start, cutoff_relative = 0.05,
                                          transcript_or_counts="NumTranscripts")
            pivoted.plot(kind='bar', stacked=True, color = sns_palette)
            locs, labels = plt.xticks()
            ax.set_xticks(ticks = locs)
            ax.set_xticklabels(labels = [label.get_text()[0:20] for label in labels])
            plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.,title=None)
            plt.tight_layout()
            os.system("mkdir -p " + results_viz_dir)
            plt.savefig(os.path.join(results_viz_dir, l + '_transcripts.png'),dpi=100)
            plt.show()
            plt.close()
        else:
            f, (ax1, ax2) = plt.subplots(1, 2, sharey=True, figsize=(15,7.5))
            ax1.set_facecolor('white')
            ax2.set_facecolor('white')
            pivoted = createPlotDataFrame(curr_df_start, cutoff_relative = 0.05,
                                          transcript_or_counts="NumTranscripts")
            pivoted.plot(kind='bar', stacked=True, width=1, color = sns_palette,
                         title="Transcripts", ax = ax1)
            locs = ax1.get_xticks()
            labels = ax1.get_xticklabels()
            ax1.set_xticks(locs)
            ax1.set_xticklabels([label.get_text()[0:20] for label in labels])
            ax1.legend(bbox_to_anchor=(1.1, 1.1),title=None)
            handles, labels = ax1.get_legend_handles_labels()
            lgd = ax1.legend(handles, labels, loc='upper center',\
                             bbox_to_anchor=(0.5,-0.1),bbox_inches='tight', title=None)
            pivoted = createPlotDataFrame(curr_dfs_start, cutoff_relative = 0.05,
                                          transcript_or_counts="Counts")
            pivoted.plot(kind='bar', stacked=True, width=1, color = sns_palette,
                         title="Counts", ax = ax2)
            locs = ax2.get_xticks()
            labels = ax2.get_xticklabels()
            ax2.set_xticks(ticks = locs)
            ax2.set_xticklabels(labels = [label.get_text()[0:20] for label in labels])
            ax2.legend(bbox_to_anchor=(1.1, 1.1),title=None)
            handles, labels = ax2.get_legend_handles_labels()
            lgd = ax2.legend(handles, labels, loc='upper center',\
                             bbox_to_anchor=(0.5,-0.1),bbox_inches='tight')
            plt.tight_layout()
            os.system("mkdir -p " + results_viz_dir)
            plt.savefig(os.path.join(results_viz_dir, l +\
                                     '_counts_and_transcripts.png'),dpi=100)
            plt.show()
            plt.close()
