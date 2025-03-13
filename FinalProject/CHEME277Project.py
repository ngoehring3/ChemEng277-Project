import os
import sys
import json
import copy
import pandas
import numpy as np
from sklearn.linear_model import LinearRegression, Ridge, RidgeCV
from sklearn.cluster import KMeans
from sklearn.preprocessing import StandardScaler

#Refine Data Structure
def data_refine(data, filename):
    hold = []
    metadict = {}
    large = data["array"]
    count = 0
    metadict['name'] = large[0][0]['substance']['name']
    metadict['study'] = data['study']['sid']
    (drug, dose) = filename.split('_', 1)
    (dose, fi) = dose.split('.', 1)
    metadict['dose'] = dose
    metadict['tissue'] = large[0][0]['tissue']['name']
    if 'group' in large[0][0]:
        metadict['count'] = large[0][0]['group']['count']
    else:
        metadict['count'] = 1
    metadict['measurment'] = large[0][0]['measurement_type']['name']
    metadict['x_unit'] = large[0][0]['time_unit']
    metadict['y_unit'] = large[0][0]['unit']
    x_value = []
    y_value = []
    sd_value = []
    se_value = []
    cv_value = []
    for x in large:
        if large[count][0]['time_unit'] == 'hr':
            x_value.append(large[count][0]['time']*60)
            metadict['x_unit'] = 'min'
        else:
            x_value.append(large[count][0]['time'])
        if 'mean' in large[count][0]:
            y_value.append(large[count][0]['mean'])
        if 'median' in large[count][0]:
            y_value.append(large[count][0]['median'])
        if 'sd' in large[count][0]:
            sd_value.append(large[count][0]['sd'])
        if 'cv' in large[count][0]:
            cv_value.append(large[count][0]['cv'])
        if 'se' in large[count][0]:
            cv_value.append(large[count][0]['se'])
        if 'value' in large[count][0]:
            y_value.append(large[count][0]['value'])

        if 'value' not in large[count][0] and 'median' not in large[count][0] and 'mean' not in large[count][0]:
            x_value.pop()
        count += 1
    metadict['x_value'] = x_value
    metadict['y_value'] = y_value
    if len(metadict['y_value']) != len(metadict['x_value']):
        print(metadict['name'])
        print(metadict['x_value'])
        print(metadict['y_value'])
        print()
    metadict['sd_value'] = sd_value
    metadict['cv_value'] = cv_value
    metadict['se_value'] = se_value
    return metadict

#Import Data
def load_json_files_from_folder(folder_path):
    end_data = {}
    drug_props = []
    for filename in os.listdir(folder_path):
        if filename.endswith('.json'):
            file_path = os.path.join(folder_path, filename)
            with open(file_path, 'r') as file:
                data = json.load(file)
                refined_data = data_refine(data, filename)
                if refined_data['name'] not in end_data:
                    end_data[refined_data['name']] = []
                if data['study']['sid'] not in end_data:
                    end_data[refined_data['name']].append(refined_data)
        if filename.endswith('.xlsx'):
            file_path = os.path.join(folder_path, filename)
            df = pandas.read_excel(file_path)
            data_list = df.to_dict(orient='records')
            drug_props.append(data_list)
    return end_data, drug_props

def time_index_total_exp(study, time):
    count = 0
    for x in study['x_value']:
        if int(x) > int(time) or count == len(study['x_value'])-1:
            return count
        count += 1

def PK_data_calculation(PK_database, time):

    calcsdict = {}
    for drug in PK_database.keys():
        if drug not in calcsdict:
            calcsdict[drug] = {}
        for study in PK_database[drug]:

            #determining total exposure
            if study['dose'] not in calcsdict[drug]:
                calcsdict[drug][study['dose']] = {}
            if 'total_exp' not in calcsdict[drug][study['dose']]:
                calcsdict[drug][study['dose']]['total_exp'] = []
            if 'max_exp' not in calcsdict[drug][study['dose']]:
                calcsdict[drug][study['dose']]['max_exp'] = []
            if 'half_life' not in calcsdict[drug][study['dose']]:
                calcsdict[drug][study['dose']]['half_life'] = []
            if study['y_unit'] == 'gram / liter':
                count = 0
            index = time_index_total_exp(study, time)
            for x in study['x_value']:
                AUC = np.trapezoid(study['y_value'][0:index - 1], study['x_value'][0:index - 1])
            if AUC != 0:
                calcsdict[drug][study['dose']]['total_exp'].append((study['count'], AUC*1000))
            else:
                calcsdict[drug][study['dose']]['total_exp'].append(None)

            #determine max exposure
            max_exp = max(study['y_value'])
            max_index = study['y_value'].index(max_exp)
            max_time = study['x_value'][max_index]
            calcsdict[drug][study['dose']]['max_exp'].append((study['count'], max_exp*1000, max_time))
            #remove outliers?

            #determine half life
            prox = []
            half_life_y = np.array([[max_exp/2]])
            post_max_y = study['y_value'][max_index:]
            post_max_x = study['x_value'][max_index:]
            x = np.array(post_max_x)
            y = np.array(post_max_y)
            x_array = []
            for z in x:
                x_array.append([z])
            y_array = []
            for t in y:
                y_array.append([t])
            if len(y_array) == len(x_array):
                mod = LinearRegression().fit(y_array, x_array)
                half_life = mod.predict(half_life_y)
                half_life = int(half_life[0][0])
                calcsdict[drug][study['dose']]['half_life'].append((study['count'], half_life))
            else:
                calcsdict[drug][study['dose']]['total_exp'].append(None)
    return calcsdict

def PK_data_consolidation(calc_data, Drug_Prop_matrix):
    total_drug_data = {}
    drug_data = {}
    for drug in calc_data:
        num_data = []
        total_num_data = []
        dose_num = 0
        total_dose_num = 0
        for dose in calc_data[drug]:
            num_dose = int(dose)
            total_data =[]
            total_data.extend(Drug_Prop_matrix[total_dose_num][1:])
            total_data.append(num_dose)
            dose_data = [num_dose]
            for data in calc_data[drug][dose]:
                con_data = 1
                n = 0
                for value in calc_data[drug][dose][data]:
                    if value != None:
                           #geometric mean to reduce outliers
                           norm_data = value[1] ** value[0]
                           con_data *= norm_data
                           n += value[0]
                if value != None and n != 0:
                    geo_mean = pow(con_data, 1 / n)
                else:
                    geo_mean = None
                calc_data[drug][dose][data] = geo_mean
                dose_data.append(geo_mean)
                total_data.append(geo_mean)
            num_data.append(dose_data)
            if None not in total_data:
                total_num_data.append(total_data)
            total_dose_num += 1
        if drug not in drug_data:
            drug_data[drug] = []
            dose_num += 1
        if drug not in total_drug_data:
            total_drug_data[drug] = []
            total_dose_num += 1
        drug_data[drug] = num_data
        total_drug_data[drug] = total_num_data
    return calc_data,drug_data,total_drug_data

def Drug_Prop_consolidation(drugs):
    #pull values from dictionary
    drug_values = []
    count = 1
    for data in drugs[0]:
        values = []
        for item in data:
            val = data[item]
            if item == 'Drug':
                val = count
            values.append(val)
        drug_values.append(values)
        count += 1
    return(drug_values)

def Drug_Prop_clustering(Drug_Prop, clusters):
    X = []
    Y = []
    x = 0
    for y in range(len(Drug_Prop)):
        X.append(Drug_Prop[y][x])
        Y.append(Drug_Prop[y][1:])
    kmeans_drugs = {'init': 'random', 'n_init': 10, 'random_state': 1}
    scaled_props = StandardScaler().fit_transform(Y)
    kmeans = KMeans(n_clusters=clusters, **kmeans_drugs)
    kmeans.fit(scaled_props)
    drug_labels = kmeans.labels_

    #dict of indexies
    indx_labels = {}
    index = 0
    for label in drug_labels:
        if label not in indx_labels:
            indx_labels[label] = []
        indx_labels[int(label)].append(index)
        index += 1
    return drug_labels, indx_labels
def mean_data_grouping(drug_labels, drug_matrix, clusters):
    #group indecies
    lst_groups = {}
    for num in range(clusters):
        for group in drug_labels:
            if group not in lst_groups:
                lst_groups[group] = []
            indexs = np.where(drug_labels == group)[0]
            lst_groups[group] = indexs
    #add data to means index and find mean of grouping
    lst_data_group = []
    for num in lst_groups:
        data = []
        indexs = lst_groups[num]
        for index in indexs:
            data.append(drug_matrix[index])
        mean_mat = []
        mean_mat.append(int(num))
        t_range = range(len(drug_matrix[0]))
        for t in t_range[1:]:
            sum_data = 0
            for z in range(len(data)):
                sum_data += data[z][t]
            mean_data = sum_data / len(data)
            mean_mat.append(mean_data)
        lst_data_group.append(mean_mat)
        lst_data_group.sort(key=lambda x:x[0])
    return(lst_data_group, lst_groups)

def distance_euclidian(test, train):
    sum_dis_total = 0
    for d in range(len(test)):
        sum_dis = pow((train[d] - test[d]), 2)
        sum_dis_total += sum_dis
    return pow(sum_dis_total, 0.5)

def Test_Data_grouping(group_means_prop, Drug_Prop_matrix_test):
    count_test = 0
    grouping = []
    for test in Drug_Prop_matrix_test:
        dist_data_test = Drug_Prop_matrix_test[count_test][1:]
        dist_total = []
        count_train = 0
        for train in group_means_prop:
            dist_data_train = group_means_prop[count_train][1:]
            euc_dist = distance_euclidian(dist_data_test, dist_data_train)
            dist_total.append(euc_dist)
            count_train += 1
        group = dist_total.index(min(dist_total))
        count_test += 1
        grouping.append(group)
    return grouping

def model_data_for_cluster(drug_index, PK_matrix, test_groups, test_matrix):
    test_train_data = []
    score_ridge = []
    test_count = 0
    for test in test_groups:
        train_index = drug_index[test]
        train_data = []
        pk_keys = list(PK_matrix.keys())
        for i in train_index:
            train_data.extend(PK_matrix[pk_keys[i]])
        X = []
        Y = []
        for y in range(len(train_data)):
            x = train_data[y][:len(train_data[0])-3]
            y = train_data[y][len(train_data[0])-3:]
            X.append(x)
            Y.append(y)
        ridge = RidgeCV()
        rid_reg = ridge.fit(X, Y)
        score_reg = ridge.score(X, Y)
        pred_drug_prop = rid_reg.predict([test_matrix[test_count][1:]])
        test_train_data.append(pred_drug_prop)
        score_ridge.append(score_reg)
        test_count += 1
    return test_train_data, score_ridge

def main():
    args = sys.argv[1:]
    if len(args) != 6:
        return print("Incorrect number of arguments")
    Train_File = args[0]
    Test_File = args[1]
    Result_Train = args[2]
    Result_Test = args[3]
    Time_cutoff = args[4]
    Cluster_num = int(args[5])

    """ Import, structure, and analyse PK data """
    #import and structure PK training data, import drug properties from train folder and builds a dict to hold them.
    (PK_database, Drug_Prop) = load_json_files_from_folder(Train_File)

    #calculate total exposure, max exposure, and half life per dose based on PK data given
    Calc_Data= PK_data_calculation(PK_database, Time_cutoff)

    #Take the geometric mean for each drug dose and create a dict with each drug dose, matrix of averaged PK data,
    # and a matrix with properties plus PK data
    Drug_Prop_matrix = Drug_Prop_consolidation(Drug_Prop)
    (Data_per_dose, PK_matrix, Total_matrix) = PK_data_consolidation(Calc_Data, Drug_Prop_matrix)

    """ Pull test data from test folder and create the same data structures for test data"""
    (PK_database_test, Drug_Prop_test) = load_json_files_from_folder(Test_File)
    Drug_Prop_matrix_test = Drug_Prop_consolidation(Drug_Prop_test)
    Test_Calc_Data = PK_data_calculation(PK_database_test, Time_cutoff)
    (Data_per_dose_test, PK_matrix_test, Total_matrix_test) = PK_data_consolidation(Test_Calc_Data, Drug_Prop_matrix_test)

    """ Manipulate the test data to etract predicted dose and reduce matrix to match training data """
    test_matrix = copy.deepcopy(Drug_Prop_matrix_test)
    #grab simulated dose
    test_dose = []
    for x in range(len(Drug_Prop_matrix_test)):
        dose = Drug_Prop_matrix_test[x].pop()
        test_dose.append(dose)

    """ K-means clustering of drug propeties"""
    #K-means clustering
    drug_labels, drug_index = Drug_Prop_clustering(Drug_Prop_matrix, Cluster_num)
    #calculate means for testing
    (group_means_prop, group_means_index) = mean_data_grouping(drug_labels, Drug_Prop_matrix, Cluster_num)

    """ Group test data into the k-means"""
    test_groups = Test_Data_grouping(group_means_prop, Drug_Prop_matrix_test)

    """ Train model for each drug based on dose in cluster; test against train data"""
    #prep of train data to use on model to match dimensions of test data
    dose_lst_train = []
    drug_labels_train = []
    train_dose = []
    drug_count = 0
    for drug in Calc_Data:
        dose = list(Calc_Data[drug].keys())
        train_dose.append(dose)
        for item in dose:
            drug_labels_train.extend([drug_labels[drug_count]])
            data = copy.deepcopy(Drug_Prop_matrix[drug_count])
            data.extend([int(item)])
            dose_lst_train.append(data)
        drug_count += 1

    #running train data predictions
    predicted_values_train, model_score_train = model_data_for_cluster(drug_index, Total_matrix, drug_labels_train, dose_lst_train)
    #running test data predictions
    predicted_values_test, model_score_test = model_data_for_cluster(drug_index, Total_matrix, test_groups, test_matrix)

    """ Append predicted values to the drug values per dose dictionary"""
    #export training data
    count = 0
    drug_count_x = 0
    train_final_key = list(Data_per_dose.keys())
    for dose_x in train_dose:
        for d in train_dose[count]:
            Data_per_dose[train_final_key[count]][str(d)]['cluster_num'] = drug_labels[count]
            Data_per_dose[train_final_key[count]][str(d)]['total_exp_pred'] = predicted_values_train[drug_count_x][0][0]
            Data_per_dose[train_final_key[count]][str(d)]['max_exp_pred'] = predicted_values_train[drug_count_x][0][1]
            Data_per_dose[train_final_key[count]][str(d)]['half_life_pred'] = predicted_values_train[drug_count_x][0][2]
            drug_count_x += 1
        count += 1
    #export training data
    count = 0
    train_final_key = list(Data_per_dose_test.keys())
    for dose_x in test_dose:
        Data_per_dose_test[train_final_key[count]][str(dose_x)]['cluster_num'] = test_groups[count]
        Data_per_dose_test[train_final_key[count]][str(dose_x)]['total_exp_pred'] = predicted_values_test[count][0][0]
        Data_per_dose_test[train_final_key[count]][str(dose_x)]['max_exp_pred'] = predicted_values_test[count][0][1]
        Data_per_dose_test[train_final_key[count]][str(dose_x)]['half_life_pred'] = predicted_values_test[count][0][2]
        count += 1

    """ Write all results to excel files"""
    #export train data
    for drug in Data_per_dose:
        file_train = pandas.DataFrame(Data_per_dose[drug])
        name = drug + '.xlsx'
        path = Result_Train
        name_path = path + "\\" + name
        file_train.to_excel(name_path, index=True)
    #export test data
    for drug in Data_per_dose_test:
        file_train = pandas.DataFrame(Data_per_dose_test[drug])
        name = drug + '.xlsx'
        path = Result_Test
        name_path = path + "\\" + name
        file_train.to_excel(name_path, index=True)
    #export cluster means
    file_train = pandas.DataFrame(group_means_prop)
    name = 'cluster_properties.xlsx'
    path = Result_Train
    name_path = path + "\\" + name
    file_train.to_excel(name_path, index=True)

if __name__ == "__main__":
    main()